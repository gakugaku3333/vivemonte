"""線源スペクトル生成 — SpekPy統合とヒール効果の軸外スペクトル。

スペクトルはSpekPyで生成する（タングステン陽極、カサレイ物理モデルの
SpekPy既定値。陽極角は scene.source.anode_angle_deg、既定12度）。
SpekPy未インストール環境ではKramers則＋Al濾過減弱の粗い近似にフォールバック
する。explicit な scene.source.spectrum（[{energy_keV, weight}, ...]）が
あればどちらより優先する。
"""
from __future__ import annotations

import warnings as _warnings

import numpy as np

from .materials import linear_mu

try:
    import spekpy as _spekpy
    _HAS_SPEKPY = True
except ImportError:
    _spekpy = None
    _HAS_SPEKPY = False


_spectrum_cache: dict = {}
_heel_cache: dict = {}


def _spekpy_spectrum(kvp: float, filtration_mm_al: float, anode_angle_deg: float):
    """SpekPyでのスペクトル生成（プロセスローカルなdictキャッシュ）。

    functools.lru_cacheではなく明示的なdictを使うのは、並列ワーカー
    起動時にSpekPy呼び出し（1回あたり約0.9秒、chest_room実測）を
    親プロセスで1回だけ行い、export_caches/import_caches経由でワーカーへ
    転送してキャッシュを事前に温めるため（docs/plan_phase3_parallel.md
    「積み残し」節）。lru_cacheは外部からのシード（値の事前挿入）が
    正式なAPIとしてサポートされていないためdictに切り替えた。
    """
    key = (kvp, filtration_mm_al, anode_angle_deg)
    if key not in _spectrum_cache:
        s = _spekpy.Spek(kvp=kvp, th=anode_angle_deg)
        s.filter("Al", filtration_mm_al)
        e_mid, phi = s.get_spectrum(edges=False)
        w = np.clip(np.asarray(phi, dtype=float), 0.0, None)
        _spectrum_cache[key] = (np.asarray(e_mid, dtype=float), w / w.sum())
    return _spectrum_cache[key]


def _kramers_fallback_spectrum(kvp: float, filtration_mm_al: float, n_bins: int = 60):
    e = np.linspace(5.0, kvp, n_bins + 1)
    e_mid = 0.5 * (e[:-1] + e[1:])
    raw = np.clip(e_mid * (kvp - e_mid), 0, None)  # Kramers則（未濾過、特性X線なし）
    mu_al = linear_mu("aluminum", e_mid)
    atten = np.exp(-mu_al * filtration_mm_al / 10.0)
    w = raw * atten
    return e_mid, w / w.sum()


def _default_spectrum(kvp: float, filtration_mm_al: float, anode_angle_deg: float = 12.0):
    if _HAS_SPEKPY:
        return _spekpy_spectrum(float(kvp), float(filtration_mm_al), float(anode_angle_deg))
    _warnings.warn("spekpy が見つからないためKramers則の粗い近似スペクトルを使用します。"
                    "`pip install spekpy` を推奨します。", stacklevel=2)
    return _kramers_fallback_spectrum(kvp, filtration_mm_al)


def sample_spectrum(src: dict, n: int, rng: np.random.Generator) -> np.ndarray:
    spec = src.get("spectrum")
    if spec:
        e = np.array([s["energy_keV"] for s in spec], dtype=float)
        w = np.array([s["weight"] for s in spec], dtype=float)
        w = w / w.sum()
    else:
        e, w = _default_spectrum(src["kvp"], src.get("filtration_mm_al", 2.5),
                                  src.get("anode_angle_deg", 12.0))
    return e[rng.choice(len(e), size=n, p=w)]


_HEEL_N_BINS = 15


def _heel_spectra(kvp: float, filtration_mm_al: float, anode_angle_deg: float,
                   sid: float, span_cm: float, n_bins: int = _HEEL_N_BINS):
    """ヒール軸ビンごとの(中心座標, スペクトル列, 相対フルエンス, 絶対フルエンス/mAs)。

    座標sはSID面上のヒール軸方向オフセット[cm]で、s>0が陽極側。
    SpekPyの軸外計算（座標系はx<0が陽極側 — 実測確認済み）で各ビンの
    スペクトルとフルエンスを求める。陽極カットオフ（take-off角が0以下に
    なる領域）ではフルエンス0になり、その旨警告する。
    """
    key = (kvp, filtration_mm_al, anode_angle_deg, sid, span_cm, n_bins)
    if key in _heel_cache:
        return _heel_cache[key]
    if not _HAS_SPEKPY:
        raise RuntimeError("ヒール効果の計算にはspekpyが必要です"
                            "（`.venv/bin/pip install spekpy`）")
    centers = (np.arange(n_bins) + 0.5) / n_bins * 2.0 * span_cm - span_cm
    spectra = []
    flu_abs = []
    for s in centers:
        sp = _spekpy.Spek(kvp=kvp, th=anode_angle_deg, z=sid, x=-s)  # s>0(陽極側)→SpekPy -x
        sp.filter("Al", filtration_mm_al)
        with _warnings.catch_warnings():
            _warnings.simplefilter("ignore")  # カットオフ域のゼロ放出警告は下でまとめて出す
            e_mid, phi = sp.get_spectrum(edges=False)
            w = np.clip(np.asarray(phi, dtype=float), 0.0, None)
            tot = float(w.sum())
            flu = float(sp.get_flu()) if tot > 0 else 0.0
        spectra.append((np.asarray(e_mid, dtype=float),
                        w / tot if tot > 0 else w))
        flu_abs.append(flu)
    flu_abs = np.array(flu_abs)
    if flu_abs.max() <= 0:
        raise ValueError("ヒール効果計算: 全ビンでX線放出がゼロです。"
                          "陽極角と照射野の組み合わせを確認してください")
    if (flu_abs <= 0).any():
        _warnings.warn(
            "照射野の陽極側端が陽極カットオフ（take-off角0度）を超えており、"
            "その領域の放出はゼロです。実機でも照射できない配置なので、"
            "照射野サイズまたは陽極角を見直してください。", stacklevel=2)
    result = (centers, tuple(spectra), flu_abs / flu_abs.max(), flu_abs)
    _heel_cache[key] = result
    return result


def export_caches() -> tuple[dict, dict]:
    """プロセスローカルなスペクトルキャッシュのスナップショットを返す。

    並列ワーカー起動前に親プロセスで一度スペクトルをウォームアップし
    （SpekPy呼び出しはkvp単発で約0.9秒、ヒール効果はビンの数だけ
    さらに重い）、この結果をワーカーへpickleで転送して`import_caches`で
    再注入するために使う（docs/plan_phase3_parallel.md「積み残し」節）。
    """
    return dict(_spectrum_cache), dict(_heel_cache)


def import_caches(spectrum_cache: dict, heel_cache: dict) -> None:
    """`export_caches`で得たキャッシュ内容を現在のプロセスへ注入する。"""
    _spectrum_cache.update(spectrum_cache)
    _heel_cache.update(heel_cache)


def heel_spectra_for_source(src: dict, sid: float, span_cm: float):
    """線源設定srcからヒール軸ビンのスペクトル群を取得する（_heel_spectra参照）。"""
    return _heel_spectra(float(src["kvp"]), float(src.get("filtration_mm_al", 2.5)),
                          float(src.get("anode_angle_deg", 12.0)), float(sid), float(span_cm))
