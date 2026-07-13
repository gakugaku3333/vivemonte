"""NIST XAAMDI（X-Ray Mass Attenuation Coefficients, Hubbell & Seltzer）を取得し
vivemonte/data/nist_xaamdi/*.csv として同梱する。

μen/ρ（質量エネルギー吸収係数）の一次ソース。カーマ・吸収線量タリーに使用。
xraylib の CS_Energy は本テーブルと最大約17%乖離することを確認済みのため、
μen/ρ は必ずこの同梱データを使う。

実行: .venv/bin/python scripts/fetch_nist_xaamdi.py
出典: https://physics.nist.gov/PhysRefData/XrayMassCoef/ (public domain)
"""
from __future__ import annotations

import re
import time
import urllib.request
from pathlib import Path

BASE = "https://physics.nist.gov/PhysRefData/XrayMassCoef"
OUT = Path(__file__).resolve().parent.parent / "vivemonte" / "data" / "nist_xaamdi"

# 短縮名 → (URL片, 説明)
TARGETS = {
    # 元素
    "z01_H": ("ElemTab/z01.html", "Hydrogen"),
    "z06_C": ("ElemTab/z06.html", "Carbon"),
    "z07_N": ("ElemTab/z07.html", "Nitrogen"),
    "z08_O": ("ElemTab/z08.html", "Oxygen"),
    "z13_Al": ("ElemTab/z13.html", "Aluminum"),
    "z20_Ca": ("ElemTab/z20.html", "Calcium"),
    "z26_Fe": ("ElemTab/z26.html", "Iron"),
    "z29_Cu": ("ElemTab/z29.html", "Copper"),
    "z74_W": ("ElemTab/z74.html", "Tungsten"),
    "z82_Pb": ("ElemTab/z82.html", "Lead"),
    # 化合物・混合物
    "air": ("ComTab/air.html", "Air, Dry (near sea level)"),
    "water": ("ComTab/water.html", "Water, Liquid"),
    "soft_tissue": ("ComTab/tissue.html", "Tissue, Soft (ICRU-44)"),
    "bone": ("ComTab/bone.html", "Bone, Cortical (ICRU-44)"),
    "lung": ("ComTab/lung.html", "Lung Tissue (ICRU-44)"),
    "muscle": ("ComTab/muscle.html", "Muscle, Skeletal (ICRU-44)"),
    "adipose": ("ComTab/adipose.html", "Adipose Tissue (ICRU-44)"),
    "pmma": ("ComTab/pmma.html", "PMMA (Lucite/Perspex)"),
    "concrete": ("ComTab/concrete.html", "Concrete, Ordinary"),
    "lead_glass": ("ComTab/glass.html", "Glass, Lead"),
}

ROW_RE = re.compile(
    r"(\d\.\d+E[+-]\d+)\s+(\d\.\d+E[+-]\d+)\s+(\d\.\d+E[+-]\d+)")


def fetch_one(name: str, url_part: str, desc: str) -> int:
    url = f"{BASE}/{url_part}"
    req = urllib.request.Request(url, headers={"User-Agent": "viveMonte-data-fetch/0.1"})
    with urllib.request.urlopen(req, timeout=30) as r:
        html = r.read().decode("ascii", errors="replace")
    rows = ROW_RE.findall(html)
    if len(rows) < 10:
        raise RuntimeError(f"{name}: テーブル行の抽出に失敗（{len(rows)}行）: {url}")
    out = OUT / f"{name}.csv"
    prev_e = -1.0
    lines = [f"# {desc} — NIST XAAMDI (Hubbell & Seltzer), source: {url}",
             "# energy_keV,mu_rho_cm2_g,muen_rho_cm2_g"]
    for e_mev, mu, muen in rows:
        e = float(e_mev) * 1000.0
        # 吸収端では同一エネルギーが2行並ぶ → 補間のため僅かにずらす
        if e <= prev_e:
            e = prev_e * (1 + 1e-6)
        prev_e = e
        lines.append(f"{e:.9g},{mu},{muen}")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, (url_part, desc) in TARGETS.items():
        n = fetch_one(name, url_part, desc)
        print(f"{name:<14} {n:>3}行  ({desc})")
        time.sleep(0.5)  # NISTサーバーへの礼儀
    print(f"\n保存先: {OUT}")


if __name__ == "__main__":
    main()
