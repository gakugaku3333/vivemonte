"""Step 0机上検算: 自由電子KN vs S(q)重み付きコンプトンの平均エネルギー移行率。

plan_compton_transfer_check.md のStep 0。EGS5(INCOH=0, IBOUND=1)は全断面積のみ
束縛補正しコンプトン散乱光子のエネルギー抽選は自由電子Klein-Nishinaのままである一方、
ChatCarloが使うNIST XAAMDIのμen/ρはS(q)補正込みの散乱理論で評価されている。
この平均エネルギー移行率tの差が、Phase 2bで観測された−1.71%系統差の主因かを
EGS5再実行なしで検算する。

実行: PYTHONPATH=. .venv/bin/python docs/egs5_crosscheck/check_compton_transfer.py
"""
from __future__ import annotations

import math

import numpy as np
import xraylib

MEC2_KEV = 510.99895
HC_KEV_A = 12.398  # keV*Angstrom (hc)

# Phase 2bで実測済みのPEGS5抽出値(60 keV水, IBOUND=1) — extract_pegs5_photo_xs.py
PEGS5_PHOTO = 0.014987   # cm2/g
PEGS5_COMPTON = 0.176906  # cm2/g

# ChatCarlo側の基準値
MUEN_XAAMDI = 0.03190  # cm2/g, NIST XAAMDI 60 keV水, 格子点そのまま(PDD_RESULTS.mdで確認済み)

# 水の組成 (H2O): 元素ごとの (Z, 原子数/分子, 原子量)
WATER_ATOMS = [(1, 2, 1.008), (8, 1, 15.999)]


def kn_diff_cross_section(theta: np.ndarray, k: float) -> np.ndarray:
    """Klein-Nishina微分断面積 dsigma/dOmega [cm^2/electron/sr]。re^2倍は後で掛ける。"""
    re2 = 2.8179403262e-13 ** 2
    ratio = 1.0 / (1.0 + k * (1.0 - np.cos(theta)))  # E'/E0
    return 0.5 * re2 * ratio ** 2 * (ratio + 1.0 / ratio - np.sin(theta) ** 2)


def q_invA(theta: np.ndarray, energy_keV: float) -> np.ndarray:
    """運動量移行変数 x=sin(theta/2)/lambda [Angstrom^-1] (xraylib SF_Compt(Z,q)の引数)。"""
    lam_A = HC_KEV_A / energy_keV
    return np.sin(theta / 2.0) / lam_A


def average_transfer_fraction(energy_keV: float, use_sq: bool, n_theta: int = 20000) -> tuple[float, float]:
    """自由電子KN(use_sq=False)またはS(q)重み付き(use_sq=True)の
    (mu_compton/rho [cm2/g], 平均エネルギー移行率 t=<T>/E0) を返す。
    水の質量分率で元素合成する。
    """
    k = energy_keV / MEC2_KEV
    theta = np.linspace(1e-6, math.pi - 1e-6, n_theta)
    sin_theta = np.sin(theta)
    kn = kn_diff_cross_section(theta, k)  # cm2/electron/sr
    e_prime_over_e0 = 1.0 / (1.0 + k * (1.0 - np.cos(theta)))
    t_theta = 1.0 - e_prime_over_e0  # T/E0 = 1 - E'/E0

    NA = 6.02214076e23
    total_A = sum(n * A for _, n, A in WATER_ATOMS)  # H2Oのモル質量相当

    sigma_sum = np.zeros_like(theta)      # dsigma/dOmega [cm2/g/sr] 積算用(水1g換算)
    sigma_t_sum = np.zeros_like(theta)    # 上記×T(theta)

    for Z, n_atoms, A in WATER_ATOMS:
        n_per_g = NA * n_atoms / total_A  # その元素の原子数密度[1/g]
        if use_sq:
            q = q_invA(theta, energy_keV)
            # xraylibのSF_Comptテーブルはq<~1e-3(元素依存)で外挿エラーになるためクリップ。
            # S(q)は物理的にq->0で0に漸近する領域なので、クリップは前方散乱の寄与
            # (元々ほぼゼロ)を歪めない。
            q_clipped = np.clip(q, 2e-3, None)
            sq = np.array([xraylib.SF_Compt(Z, float(qi)) for qi in q_clipped])
            sq = np.where(q < 2e-3, 0.0, sq)
            weight = sq  # 電子1個あたりのKNにS(q,Z)を掛けて原子あたりに変換(Z電子分をS(q)が代表)
        else:
            weight = float(Z)  # 束縛効果なし: Z個の自由電子として単純加算
        diff = n_per_g * weight * kn  # [cm2/g/sr]
        sigma_sum += diff
        sigma_t_sum += diff * t_theta

    # 立体角積分: 2*pi*sin(theta) dtheta
    mu_rho = 2 * math.pi * np.trapezoid(sigma_sum * sin_theta, theta)
    mu_rho_t = 2 * math.pi * np.trapezoid(sigma_t_sum * sin_theta, theta)
    t_avg = mu_rho_t / mu_rho
    return mu_rho, t_avg


def main() -> None:
    E = 60.0

    print(f"### Step 0: 60 keV水、コンプトン平均エネルギー移行率の三者比較 ###\n")

    mu_kn, t_kn = average_transfer_fraction(E, use_sq=False)
    mu_sq, t_sq = average_transfer_fraction(E, use_sq=True)

    photo_xr = xraylib.CS_Photo_CP("Water, Liquid", E)
    compt_xr = xraylib.CS_Compt_CP("Water, Liquid", E)
    compton_transfer_xaamdi = MUEN_XAAMDI - photo_xr
    t_xaamdi = compton_transfer_xaamdi / compt_xr

    print(f"自由電子KN(束縛なし)      : mu_compton/rho={mu_kn:.6f} cm2/g  t_KN={t_kn:.6f}")
    print(f"S(q)重み付き(束縛あり)     : mu_compton/rho={mu_sq:.6f} cm2/g  t_Sq={t_sq:.6f}")
    print(f"xraylib(EPDL) mu_compton/rho = {compt_xr:.6f} cm2/g  (参考: KN/Sq計算の断面積側チェック)")
    print(f"XAAMDI暗黙値              : t_XAAMDI = {t_xaamdi:.6f}  "
          f"(muen={MUEN_XAAMDI}, photo_xr={photo_xr:.6f}, compt_xr={compt_xr:.6f})")
    print()
    print(f"t_KN / t_XAAMDI = {t_kn/t_xaamdi:.4f}  ({(t_kn/t_xaamdi-1)*100:+.2f}%)")
    print(f"t_Sq / t_XAAMDI = {t_sq/t_xaamdi:.4f}  ({(t_sq/t_xaamdi-1)*100:+.2f}%)  <- 整合性チェック(0%に近いはず)")
    print()

    print("### Delta_pred: PEGS5断面積 + 各tでのkerma係数 vs ChatCarlo(XAAMDI) ###")
    for label, t in [("t_KN(EGS5相当, INCOH=0)", t_kn), ("t_Sq(束縛補正込み)", t_sq), ("t_XAAMDI(参考)", t_xaamdi)]:
        k_egs5_eff = PEGS5_PHOTO + PEGS5_COMPTON * t
        delta = (k_egs5_eff - MUEN_XAAMDI) / MUEN_XAAMDI * 100
        print(f"  {label:30s}: k_eff={k_egs5_eff:.6f} cm2/g  Delta_pred={delta:+.3f}%")

    print()
    print("観測されたPhase 2bの系統差(47ビン平均): -1.71% (範囲 -0.71%〜-3.45%)")
    print("事前登録基準: 主因=-1.2%〜-2.2%, 棄却=|Delta_pred|<0.5%, それ以外=部分的要因")

    print()
    print("### 副検算: t_KN不足のエネルギー依存性(20-60 keV) ###")
    print("(vive-auditor所見7: このスキャンをスクリプト内に残し単体再現可能にする)")
    for e in (20, 30, 40, 50, 60):
        mu_kn_e, t_kn_e = average_transfer_fraction(float(e), use_sq=False)
        mu_sq_e, t_sq_e = average_transfer_fraction(float(e), use_sq=True)
        photo_e = xraylib.CS_Photo_CP("Water, Liquid", float(e))
        compt_e = xraylib.CS_Compt_CP("Water, Liquid", float(e))
        print(f"  E={e:3d}keV  t_KN={t_kn_e:.5f}  t_Sq={t_sq_e:.5f}  "
              f"(t_KN/t_Sq-1)={(t_kn_e/t_sq_e-1)*100:+.2f}%  "
              f"photo/(photo+compton)={photo_e/(photo_e+compt_e)*100:.1f}%")


if __name__ == "__main__":
    main()
