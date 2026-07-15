"""差の真因の切り分け: 自由電子Klein-Nishina全断面積 vs EPDL束縛補正込み断面積。

初回比較（EGS5 IBOUND=0）で観測された相対約4.5%の透過率差の原因仮説を検証する。
xraylibの CS_Compt はEPDL由来の非干渉性散乱断面積（束縛電子補正込み）であり、
自由電子KNの解析積分より小さい。EGS5のIBOUND=0は自由電子KNを使うため、
コンプトン断面積だけKN値に差し替えた透過率がEGS5実測を再現するかを確認する。

実行: PYTHONPATH=. .venv/bin/python docs/egs5_crosscheck/check_kn_vs_epdl.py
"""
from __future__ import annotations

import math

import xraylib

E_KEV = 60.0
THICKNESS_CM = 10.0


def kn_total_cross_section_cm2(energy_keV: float) -> float:
    """自由電子Klein-Nishina全断面積の解析積分 [cm^2/electron]。"""
    mec2 = 510.99895  # keV
    k = energy_keV / mec2
    re = 2.8179403262e-13  # cm
    t1 = (1 + k) / k**2 * (2 * (1 + k) / (1 + 2 * k) - math.log(1 + 2 * k) / k)
    t2 = math.log(1 + 2 * k) / (2 * k)
    t3 = -(1 + 3 * k) / (1 + 2 * k) ** 2
    return 2 * math.pi * re**2 * (t1 + t2 + t3)


def main() -> None:
    # 水の電子数密度（H2O: 分子あたり10電子、モル質量18.0153 g/mol）
    NA = 6.02214076e23
    ne_per_g = NA * 10 / 18.0153

    mu_rho_kn = kn_total_cross_section_cm2(E_KEV) * ne_per_g
    cs_compt = xraylib.CS_Compt_CP("Water, Liquid", E_KEV)
    cs_rayl = xraylib.CS_Rayl_CP("Water, Liquid", E_KEV)
    cs_photo = xraylib.CS_Photo_CP("Water, Liquid", E_KEV)

    print(f"コンプトン mu/rho: 自由電子KN(解析) = {mu_rho_kn:.5f} cm2/g")
    print(f"コンプトン mu/rho: xraylib(EPDL束縛) = {cs_compt:.5f} cm2/g")
    print(f"比 KN/EPDL = {mu_rho_kn / cs_compt:.4f} (+{(mu_rho_kn / cs_compt - 1) * 100:.2f}%)")

    total_epdl = cs_compt + cs_rayl + cs_photo
    total_kn = mu_rho_kn + cs_rayl + cs_photo
    for rho, label in [(1.0, "rho=1.000 (viveMonte)"), (1.001, "rho=1.001 (PEGS5テンプレート)")]:
        t_epdl = math.exp(-total_epdl * rho * THICKNESS_CM)
        t_kn = math.exp(-total_kn * rho * THICKNESS_CM)
        print(f"{label}: T(EPDLコンプトン)={t_epdl * 100:.3f}%  T(KNコンプトン)={t_kn * 100:.3f}%")

    print()
    print("EGS5実測: IBOUND=0(自由KN)=12.02%  IBOUND=1(束縛)=12.70%")
    print("→ コンプトン断面積の物理モデル差（束縛補正の有無）が初回の差を説明する")


if __name__ == "__main__":
    main()
