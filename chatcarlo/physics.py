"""相互作用の物理サンプリング — コンプトン/レイリー/蛍光X線の角度・エネルギー抽選。

相互作用種別は光電/コンプトン/レイリーの3種（診断領域で対生成は無視）。
- 光電: 光電吸収元素をZ別光電断面積で抽選し、K殻蛍光X線を抽選する
  （`sample_fluorescence`）。放出されればエネルギーの一部を持ち出す光子として
  輸送を続け、放出されなければ従来どおりエネルギー全量をその場で局所吸収
  （電子飛程を無視するカーマ近似。README/[[lessons_learned]]参照）
  — いずれも輸送カーネル側（transport.py）で処理
- コンプトン: 束縛コンプトン散乱。元素をZ別コンプトン断面積で抽選し、Klein-Nishina
  微分断面積からKahn型棄却法で提案したε=E'/E・散乱角を、非干渉性散乱関数
  S(Z,q)/Z（xraylib.SF_Compt, EPDLベース）で追加棄却する（`sample_compton_bound`）。
  自由電子版（束縛効果なし）は`sample_klein_nishina`として比較用に残している
- レイリー: 弾性散乱、エネルギー変化なし。角度分布は原子形状因子F(Z,q)込みの
  微分断面積 dσ/dΩ ∝ (1+cos²θ)·F(Z,q)² から棄却法で抽出する（xraylib.FF_Rayl,
  EPDLベース）。化合物・混合物では質量分率×元素別レイリー断面積で
  構成元素をまず抽選してからその元素のF(Z,q)を使う
"""
from __future__ import annotations

import numpy as np
import xraylib

from .materials import (compton_element_weights, fluorescence_k_data,
                         incoherent_sq_table, material_groups,
                         photo_element_weights, rayleigh_element_weights,
                         rayleigh_form_factor_table)

_MEC2_KEV = 511.0
_HC_KEV_ANGSTROM = 12.3984193  # xraylib.MomentTransfと同じ定数（hc）
# K蛍光カットオフ: これ未満の線エネルギーは自material内mfpがμmオーダーで
# 実質局所吸収されるため、蛍光光子として生成しない（docs/plan_fluorescence.md参照）
_FLUOR_CUTOFF_KEV = 5.0


def _propose_free_electron_kn(a_p: np.ndarray, emin_p: np.ndarray, m_p: np.ndarray,
                               rng: np.random.Generator):
    """Kahn型棄却法の1トライアル分。ε=E'/E ∈ [1/(1+2α),1] をKN微分断面積に
    従って提案し、自由電子KNとしての受理可否を返す(eps_p, cos_p, accept)。

    g(ε) = 1/ε + ε - sin²θ(ε) は ε=ε_min（後方散乱）で最大値
    M = 1/ε_min + ε_min を取るため、それを一様提案の包絡線に使う。
    """
    xi1 = rng.random(len(a_p))
    xi2 = rng.random(len(a_p))
    eps_p = emin_p + xi1 * (1.0 - emin_p)
    cos_p = 1.0 - (1.0 / eps_p - 1.0) / a_p
    sin2_p = 1.0 - cos_p ** 2
    g = 1.0 / eps_p + eps_p - sin2_p
    accept = xi2 * m_p <= g
    return eps_p, cos_p, accept


def sample_klein_nishina(e_keV: np.ndarray, rng: np.random.Generator):
    """自由電子Klein-Nishina。束縛効果なし（`sample_compton_bound`参照）。"""
    alpha = e_keV / _MEC2_KEV
    eps_min = 1.0 / (1.0 + 2.0 * alpha)
    envelope = 1.0 / eps_min + eps_min
    n = len(e_keV)
    eps = np.empty(n)
    cos_theta = np.empty(n)
    pending = np.arange(n)
    while len(pending) > 0:
        eps_p, cos_p, accept = _propose_free_electron_kn(
            alpha[pending], eps_min[pending], envelope[pending], rng)
        acc = pending[accept]
        eps[acc] = eps_p[accept]
        cos_theta[acc] = cos_p[accept]
        pending = pending[~accept]
    return eps, cos_theta


def sample_compton_element(materials: np.ndarray, energies: np.ndarray,
                            rng: np.random.Generator) -> np.ndarray:
    """化合物・混合物の中で、コンプトン相互作用がどの構成元素で起きたかを抽選する。

    質量分率×元素別コンプトン断面積（束縛効果込み）で規格化した重みに従う
    （`sample_rayleigh_element`と対称的な設計、materials.py参照）。
    """
    z_chosen = np.empty(len(materials), dtype=int)
    for name, m in material_groups(materials):
        zs, w = compton_element_weights(name, energies[m])
        cumw = np.cumsum(w, axis=0)
        r = rng.random(int(np.sum(m)))
        idx = np.clip(np.sum(r[None, :] > cumw, axis=0), 0, len(zs) - 1)
        z_chosen[m] = zs[idx]
    return z_chosen


def sample_compton_bound(materials: np.ndarray, e_keV: np.ndarray,
                          rng: np.random.Generator):
    """S(Z,q)非干渉性散乱関数込みの束縛コンプトン散乱。ε=E'/Eと散乱角cosθを返す。

    元素をZ別コンプトン断面積で抽選した上で、自由電子KN(Kahn型棄却法)で
    (ε, cosθ)を提案し、S(Z,q)/Z（q=E·sin(θ/2)/hc、`incoherent_sq_table`参照）
    を追加の受理確率として棄却する。S(Z,q)はq→∞でZに単調収束するため
    0≤S(Z,q)/Z≤1が常に成り立ち、追加の包絡線調整は不要（2段の独立な棄却法の
    合成: 自由電子KN受理とS(q)/Z受理を両方満たした試行だけを採用する）。
    S(Z,q)はq→0で0に漸近するため、小角散乱（小さいq）ほど強く抑制される
    （Rayleigh散乱のF(Z,q)がq=0で最大になるのと対照的な挙動）。
    """
    z_array = sample_compton_element(materials, e_keV, rng)
    alpha = e_keV / _MEC2_KEV
    eps_min = 1.0 / (1.0 + 2.0 * alpha)
    envelope = 1.0 / eps_min + eps_min
    n = len(e_keV)
    eps = np.empty(n)
    cos_theta = np.empty(n)
    pending = np.arange(n)
    while len(pending) > 0:
        eps_p, cos_p, accept_kn = _propose_free_electron_kn(
            alpha[pending], eps_min[pending], envelope[pending], rng)

        e_p = e_keV[pending]
        z_p = z_array[pending]
        theta_p = np.arccos(np.clip(cos_p, -1.0, 1.0))
        q_p = e_p * np.sin(theta_p / 2.0) / _HC_KEV_ANGSTROM
        s_over_z = np.empty(len(pending))
        for z in set(z_p.tolist()):
            mz = z_p == z
            q_grid, s_grid = incoherent_sq_table(int(z))
            s_over_z[mz] = np.interp(q_p[mz], q_grid, s_grid) / z
        xi3 = rng.random(len(pending))
        accept = accept_kn & (xi3 <= s_over_z)

        acc = pending[accept]
        eps[acc] = eps_p[accept]
        cos_theta[acc] = cos_p[accept]
        pending = pending[~accept]
    return eps, cos_theta


def sample_photo_element(materials: np.ndarray, energies: np.ndarray,
                          rng: np.random.Generator) -> np.ndarray:
    """化合物・混合物の中で、光電相互作用がどの構成元素で起きたかを抽選する。

    `sample_compton_element`/`sample_rayleigh_element`と同型。質量分率×
    元素別光電断面積で規格化した重みに従う（蛍光X線サンプリング用）。
    """
    z_chosen = np.empty(len(materials), dtype=int)
    for name, m in material_groups(materials):
        zs, w = photo_element_weights(name, energies[m])
        cumw = np.cumsum(w, axis=0)
        r = rng.random(int(np.sum(m)))
        idx = np.clip(np.sum(r[None, :] > cumw, axis=0), 0, len(zs) - 1)
        z_chosen[m] = zs[idx]
    return z_chosen


def sample_fluorescence(materials: np.ndarray, e_keV: np.ndarray,
                         rng: np.random.Generator):
    """光電吸収イベント群に対しK殻蛍光X線の放出を抽選する。

    戻り値: (emit: bool配列, e_line: float配列[keV])。emit=Trueの光子は
    e_lineのK蛍光光子を等方放出する（呼び出し側でエネルギー/方向を書き換えて
    輸送続行）。emit=Falseは従来どおり全量その場で局所吸収する。

    手順（元素Zグループごと）:
    1. 光電断面積の元素分岐で吸収元素Zを抽選（`sample_photo_element`）
    2. E<=K吸収端 の光子はK殻蛍光を出せないため以降の判定をスキップ
    3. K殻イオン化確率 CS_Photo_Partial(Z,K,E)/CS_Photo(Z,E) で棄却
    4. K蛍光収率ω_Kで棄却
    5. 有効な4線（KL2/KL3/KM2/KM3）から発生確率で線を抽選し、
       線エネルギーが_FLUOR_CUTOFF_KEV未満なら放出しない（局所吸収扱い）
    """
    z_array = sample_photo_element(materials, e_keV, rng)
    n = len(e_keV)
    emit = np.zeros(n, dtype=bool)
    e_line = np.zeros(n)
    for z in set(z_array.tolist()):
        mz = np.where(z_array == z)[0]
        edge_keV, omega_k, line_energies, line_probs = fluorescence_k_data(int(z))
        if omega_k <= 0 or line_energies.size == 0:
            continue
        if line_energies.max() < _FLUOR_CUTOFF_KEV:
            # 軽元素はK線が全てカットオフ未満で決して放出されない。
            # xraylib.CS_Photo_Partialは軽元素・高エネルギーでスプライン
            # 外挿エラーを起こすことがあるため、その呼び出し自体を避ける。
            continue
        e_z = e_keV[mz]
        above_edge = e_z > edge_keV
        if not np.any(above_edge):
            continue
        idx = mz[above_edge]
        e_sub = e_z[above_edge]

        k_frac = np.array([
            xraylib.CS_Photo_Partial(int(z), xraylib.K_SHELL, float(e))
            / xraylib.CS_Photo(int(z), float(e))
            for e in e_sub
        ])
        is_k = rng.random(len(idx)) < k_frac
        is_radiative = rng.random(len(idx)) < omega_k
        candidate = idx[is_k & is_radiative]
        if len(candidate) == 0:
            continue

        cumw = np.cumsum(line_probs)
        r = rng.random(len(candidate))
        line_idx = np.clip(np.sum(r[None, :] > cumw[:, None], axis=0), 0, len(line_energies) - 1)
        chosen_e = line_energies[line_idx]
        above_cutoff = chosen_e >= _FLUOR_CUTOFF_KEV
        emit_idx = candidate[above_cutoff]
        emit[emit_idx] = True
        e_line[emit_idx] = chosen_e[above_cutoff]
    return emit, e_line


def isotropic_direction(n: int, rng: np.random.Generator) -> np.ndarray:
    """一様等方な単位方向ベクトルをn個抽選する（蛍光X線の放出方向用）。"""
    cos_theta = rng.uniform(-1.0, 1.0, n)
    sin_theta = np.sqrt(np.clip(1.0 - cos_theta ** 2, 0.0, None))
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    return np.column_stack([sin_theta * np.cos(phi), sin_theta * np.sin(phi), cos_theta])


def sample_rayleigh_element(materials: np.ndarray, energies: np.ndarray,
                             rng: np.random.Generator) -> np.ndarray:
    """化合物・混合物の中で、レイリー相互作用がどの構成元素で起きたかを抽選する。

    質量分率×元素別レイリー断面積で規格化した重みに従う（materials.py参照）。
    """
    z_chosen = np.empty(len(materials), dtype=int)
    for name, m in material_groups(materials):
        zs, w = rayleigh_element_weights(name, energies[m])  # w: (n_elem, sum(m))
        cumw = np.cumsum(w, axis=0)
        r = rng.random(int(np.sum(m)))
        idx = np.clip(np.sum(r[None, :] > cumw, axis=0), 0, len(zs) - 1)
        z_chosen[m] = zs[idx]
    return z_chosen


def sample_rayleigh_cos_theta(z_array: np.ndarray, e_array: np.ndarray,
                               rng: np.random.Generator) -> np.ndarray:
    """原子形状因子込みの微分断面積 (1+cos²θ)·F(Z,q)² を棄却法で抽出。

    q = E·sin(θ/2)/hc [Å⁻¹]（xraylib.MomentTransfと同じ定義）。F(Z,q)はq=0で
    Zを取り単調減少するため、g(cosθ)=(1+cos²θ)F(Z,q)²の最大値は前方散乱
    (θ=0, q=0)での 2Z² となり、これを棄却法の包絡線に使う。
    """
    n = len(z_array)
    cos_theta = np.empty(n)
    pending = np.arange(n)
    while len(pending) > 0:
        zp = z_array[pending]
        ep = e_array[pending]
        c = rng.uniform(-1.0, 1.0, len(pending))
        theta = np.arccos(c)
        q = ep * np.sin(theta / 2.0) / _HC_KEV_ANGSTROM

        f = np.empty(len(pending))
        for z in set(zp.tolist()):
            m = zp == z
            q_grid, f_grid = rayleigh_form_factor_table(int(z))
            f[m] = np.interp(q[m], q_grid, f_grid)

        g = (1.0 + c ** 2) * f ** 2
        envelope = 2.0 * zp.astype(float) ** 2
        xi2 = rng.random(len(pending))
        accept = xi2 * envelope <= g
        acc = pending[accept]
        cos_theta[acc] = c[accept]
        pending = pending[~accept]
    return cos_theta


def scatter_direction(d: np.ndarray, cos_theta: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """散乱角cosθと一様方位角から、現在の飛行方向dに対する新しい方向を返す。"""
    n = d.shape[0]
    sin_theta = np.sqrt(np.clip(1.0 - cos_theta ** 2, 0.0, None))
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    up = np.where((np.abs(d[:, 2]) < 0.999)[:, None],
                  np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0]))
    u = np.cross(up, d)
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    v = np.cross(d, u)
    new_dir = ((sin_theta * np.cos(phi))[:, None] * u
               + (sin_theta * np.sin(phi))[:, None] * v
               + cos_theta[:, None] * d)
    return new_dir / np.linalg.norm(new_dir, axis=1, keepdims=True)
