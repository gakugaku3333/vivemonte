"""線源サンプリング — 点線源＋照射野（rect/cone）、CTガントリー回転、ヒール効果。

点線源から照射野への発散ビームで光子(位置, 方向, エネルギー)を生成する。
スペクトル自体の生成は vivemonte/spectrum.py に分離してある。
mAs指定時の絶対光子数校正（photon_count_through_field）もここに置く
（照射野の幾何と面平均フルエンスの計算を共有するため）。
"""
from __future__ import annotations

import math

import numpy as np

from .spectrum import _HAS_SPEKPY, _spekpy, heel_spectra_for_source, sample_spectrum

_ROTATION_AXES = {"x": np.array([1.0, 0.0, 0.0]), "y": np.array([0.0, 1.0, 0.0]),
                   "z": np.array([0.0, 0.0, 1.0])}


def _rotate_batch(v: np.ndarray, axis: np.ndarray, angles: np.ndarray) -> np.ndarray:
    """固定ベクトルvを、軸axis周りに角度angles[i]だけ回転した結果をi行目に返す（Rodrigues）。"""
    cos_a = np.cos(angles)[:, None]
    sin_a = np.sin(angles)[:, None]
    return (v[None, :] * cos_a + np.cross(axis, v)[None, :] * sin_a
            + axis[None, :] * (axis @ v) * (1.0 - cos_a))


def beam_basis(d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """ビーム方向dに直交する局所基底(u, v)。SID面上の照射野座標に使う。"""
    if abs(d[2]) < 0.999:
        u = np.array([-d[1], d[0], 0.0])
    else:
        u = np.array([1.0, 0.0, 0.0])
    u = u / np.linalg.norm(u)
    return u, np.cross(d, u)


def cone_half_angle_rad(fld: dict) -> float:
    """cone照射野の半頂角[rad]。開口はSID面での直径diameter_cmで指定する。"""
    return math.atan((fld["diameter_cm"] / 2.0) / fld["sid_cm"])


def _heel_axis_coeffs(src: dict, d: np.ndarray, u: np.ndarray, v: np.ndarray):
    """anode_direction（世界座標）をビーム直交面に射影し、(u,v)基底の係数で返す。

    係数(h_u, h_v)は局所基底に対して定義されるため、rotation使用時も
    管球と一緒に回るヒール軸（実機の陽極-陰極軸は管球に固定）を正しく表す。
    """
    a = np.asarray(src["anode_direction"], dtype=float)
    a = a - (a @ d) * d
    norm = np.linalg.norm(a)
    if norm < 1e-9:
        raise ValueError("source.anode_direction がビーム方向と平行です"
                          "（陽極-陰極軸はビーム中心軸に直交する成分が必要）")
    a = a / norm
    return float(a @ u), float(a @ v)


def _heel_span_cm(fld: dict, h_u: float, h_v: float) -> float:
    """SID面上でのヒール軸方向の照射野半幅[cm]（軸外スペクトルのビン化範囲）。"""
    if fld.get("shape", "rect") == "cone":
        return fld["diameter_cm"] / 2.0
    w, h = fld["size_cm"]
    return abs(h_u) * w / 2.0 + abs(h_v) * h / 2.0


def sample_source_photons(src: dict, n: int, rng: np.random.Generator):
    """点線源から照射野への発散ビームで光子(位置,方向,エネルギー)を生成。

    照射野は2種類（source.field.shape）:
    - rect（既定）: SID面の矩形上に一様に目標点を抽選（コリメータ開口の
      幾何を優先する近似。立体角あたりの光子数は視野端でcos³θ分だけ過剰に
      なるが、一般撮影のSID・視野では数%）
    - cone: 中心軸周りの円錐（半頂角=atan(diameter_cm/2/sid_cm)）内で
      方向を立体角一様に抽選。実管球の等方放出に対応し、
      CBCT・広角視野で物理的により正しい

    source.heel_effect: true（+ anode_direction）でヒール効果を適用する:
    SpekPyの軸外スペクトルをヒール軸に沿ってビン化し、陽極側ほど光子数が
    少なく（棄却サンプリング）・スペクトルが硬くなる（ビン別エネルギー抽選）。
    rect/cone・rotation・ヘリカルすべてと合成可。anode_direction は
    局所基底の係数に変換されるため、rotation時は管球と一緒に回る。

    source.rotation が指定されている場合は、CTガントリー回転を光子ごとの
    角度抽選で表現する: 焦点位置・ビーム方向をisocenter周りに抽選角だけ回転させる。
    角度は既定で連続一様（実機CTの連続曝射に対応）、n_angles指定時のみ
    離散一様（トモシンセシス等の離散投影モダリティ用）。

    rotation.scan_length_cm 指定時はヘリカル撮影の位相平均近似:
    回転角と独立に、回転軸方向の一様分布で焦点を平行移動する。
    螺旋の開始位相が体に対してランダムであることを位相平均すると
    角度とテーブル位置は厳密に独立になるため、これは「位相平均された
    ヘリカル軌道」の統計的に厳密な表現である（スキャン端の半回転以内と
    over-rangingを除く）。
    """
    pos = np.asarray(src["position"], dtype=float)
    d = np.asarray(src["direction"], dtype=float)
    fld = src["field"]
    sid = fld["sid_cm"]
    u, v = beam_basis(d)

    # 焦点位置と局所基底（d,u,v）を決める。回転時は光子ごとの(n,3)配列、
    # 静止時は(3,)のままにしてブロードキャストで共通処理する。
    rot = src.get("rotation")
    if rot is not None:
        iso = np.asarray(rot["isocenter"], dtype=float)
        axis = _ROTATION_AXES[rot.get("axis", "z")]
        n_angles = rot.get("n_angles")
        if n_angles:
            angles = 2.0 * np.pi * rng.integers(0, int(n_angles), size=n) / int(n_angles)
        else:
            angles = rng.uniform(0.0, 2.0 * np.pi, size=n)
        pos_a = iso[None, :] + _rotate_batch(pos - iso, axis, angles)
        d_a = _rotate_batch(d, axis, angles)
        u_a = _rotate_batch(u, axis, angles)
        v_a = _rotate_batch(v, axis, angles)
        scan = float(rot.get("scan_length_cm") or 0.0)
        if scan > 0.0:
            shift = rng.uniform(-scan / 2.0, scan / 2.0, size=n)
            pos_a = pos_a + axis[None, :] * shift[:, None]
        origins = pos_a
    else:
        pos_a, d_a, u_a, v_a = pos, d, u, v
        origins = np.tile(pos, (n, 1))

    shape = fld.get("shape", "rect")

    # ヒール効果: 陽極-陰極軸に沿った強度低下（棄却サンプリング）と
    # 線質硬化（ビンごとのスペクトルからエネルギー抽選）を適用する
    heel = bool(src.get("heel_effect"))
    if heel:
        h_u, h_v = _heel_axis_coeffs(src, d, u, v)
        span = _heel_span_cm(fld, h_u, h_v)
        centers, spectra, rel_flu, _ = heel_spectra_for_source(src, float(sid), span)

        def _accept(draw):
            """draw(m) -> (ヒール座標s, ペイロード列...) を棄却法でn個集める。"""
            kept = None
            got = 0
            while got < n:
                m = max(int((n - got) * 1.8), 1024)
                s, *payload = draw(m)
                p = np.interp(s, centers, rel_flu)
                keep = rng.random(m) < p
                cols = [s[keep]] + [c[keep] for c in payload]
                kept = cols if kept is None else [np.concatenate([a, b])
                                                   for a, b in zip(kept, cols)]
                got = len(kept[0])
            return [c[:n] for c in kept]

    if shape == "cone":
        cos_half = np.cos(cone_half_angle_rad(fld))

        def _draw_cone(m):
            # 円錐キャップ内の立体角一様抽選: cosθ ~ U[cos(半頂角), 1]、方位角一様
            c = rng.uniform(cos_half, 1.0, m)
            ph = rng.uniform(0.0, 2.0 * np.pi, m)
            tan_t = np.sqrt(1.0 - c ** 2) / c
            s = sid * tan_t * (np.cos(ph) * h_u + np.sin(ph) * h_v) if heel else None
            return s, c, ph

        if heel:
            s_heel, cos_t, phi = _accept(_draw_cone)
        else:
            _, cos_t, phi = _draw_cone(n)
        sin_t = np.sqrt(1.0 - cos_t ** 2)
        dirs = (cos_t[:, None] * d_a
                + (sin_t * np.cos(phi))[:, None] * u_a
                + (sin_t * np.sin(phi))[:, None] * v_a)
    else:
        w, h = fld["size_cm"]

        def _draw_rect(m):
            a = rng.uniform(-w / 2, w / 2, m)
            b = rng.uniform(-h / 2, h / 2, m)
            s = a * h_u + b * h_v if heel else None
            return s, a, b

        if heel:
            s_heel, su, sv = _accept(_draw_rect)
        else:
            _, su, sv = _draw_rect(n)
        target = pos_a + d_a * sid + su[:, None] * u_a + sv[:, None] * v_a
        dirs = target - origins
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)

    if heel:
        # エネルギーは光子が属すヒールビンのスペクトルから抽選（陽極側ほど硬い）
        bin_w = 2.0 * span / len(centers)
        bins = np.clip(((s_heel + span) / bin_w).astype(int), 0, len(centers) - 1)
        # 陽極カットオフ境界の補間で受理された光子がゼロ放出ビンに落ちることが
        # あるため、非ゼロフルエンスのビン範囲にクランプする
        nz = np.where(rel_flu > 0)[0]
        bins = np.clip(bins, nz.min(), nz.max())
        energies = np.empty(n)
        for b in np.unique(bins):
            m = bins == b
            e_b, w_b = spectra[b]
            energies[m] = e_b[rng.choice(len(e_b), size=int(m.sum()), p=w_b)]
    else:
        energies = sample_spectrum(src, n, rng)
    return origins, dirs, energies


def photon_count_through_field(src: dict) -> float:
    """指定されたmAsで実際に照射野を通過する光子数（フルエンス×照射野面積）。

    SpekPyの絶対フルエンス計算（get_flu、focus-to-detector距離z=SIDでの
    photons/cm²）を使う。フィールド内のフルエンス分布は中心軸上の値で
    一様と近似する（斜入射による濾過路長の増加等は無視、教育・研究用の
    一次近似）。SpekPy未インストール環境では校正できない
    （Kramers近似フォールバックは絶対規格化を持たないため）。
    """
    if not _HAS_SPEKPY:
        raise RuntimeError("光子数校正にはspekpyが必要です（`.venv/bin/pip install spekpy`）")
    mas = src.get("mas")
    if mas is None:
        raise ValueError("source.mas が指定されていません（光子数校正には管電流時間積[mAs]が必要）")
    fld = src["field"]
    sid = fld["sid_cm"]
    if fld.get("shape", "rect") == "cone":
        area = np.pi * (fld["diameter_cm"] / 2.0) ** 2
    else:
        w, h = fld["size_cm"]
        area = w * h

    if src.get("heel_effect"):
        # ヒール適用時は中心軸値ではなく照射野の面平均フルエンスを使う
        d = np.asarray(src["direction"], dtype=float)
        u, v = beam_basis(d)
        h_u, h_v = _heel_axis_coeffs(src, d, u, v)
        span = _heel_span_cm(fld, h_u, h_v)
        if fld.get("shape", "rect") == "cone":
            g = np.linspace(-span, span, 101)
            gu, gv = np.meshgrid(g, g)
            inside = gu ** 2 + gv ** 2 <= span ** 2
            s_grid = (gu * h_u + gv * h_v)[inside]
        else:
            w0, h0 = fld["size_cm"]
            gu, gv = np.meshgrid(np.linspace(-w0 / 2, w0 / 2, 101),
                                  np.linspace(-h0 / 2, h0 / 2, 101))
            s_grid = gu * h_u + gv * h_v
        centers, _, _, flu_abs = heel_spectra_for_source(src, float(sid), span)
        flu_mean_per_mas = float(np.interp(s_grid.ravel(), centers, flu_abs).mean())
        return flu_mean_per_mas * mas * area

    s = _spekpy.Spek(kvp=src["kvp"], th=src.get("anode_angle_deg", 12.0), z=sid, mas=mas)
    s.filter("Al", src.get("filtration_mm_al", 2.5))
    fluence_per_cm2 = s.get_flu()
    return fluence_per_cm2 * area
