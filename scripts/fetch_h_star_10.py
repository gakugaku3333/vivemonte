"""光子の周辺線量当量H*(10)換算係数（ICRP Publication 74 / ICRU Report 57）を取得し
chatcarlo/data/h_star_10/photons_icrp74.csv として同梱する。

原資料（ICRP74本文）を直接パースする代わりに、OpenMCプロジェクト
（MITライセンス）が同梱するテーブル `openmc/data/dose/icrp74/photons_H10.txt`
を取得する。OpenMC側もICRP74公表値をそのまま転記したものであり、
値の由来はICRP74で変わらない。

実行: .venv/bin/python scripts/fetch_h_star_10.py
出典: https://github.com/openmc-dev/openmc/blob/develop/openmc/data/dose/icrp74/photons_H10.txt
"""
from __future__ import annotations

import re
import urllib.request
from pathlib import Path

URL = ("https://raw.githubusercontent.com/openmc-dev/openmc/develop/"
       "openmc/data/dose/icrp74/photons_H10.txt")
OUT = Path(__file__).resolve().parent.parent / "chatcarlo" / "data" / "h_star_10" / "photons_icrp74.csv"


def main() -> None:
    with urllib.request.urlopen(URL, timeout=30) as resp:
        text = resp.read().decode("utf-8")

    rows = []
    for line in text.splitlines():
        m = re.match(r"^\s*([\d.]+)\s+([\d.]+)\s*$", line)
        if m:
            e_mev, h_pSv_cm2 = float(m.group(1)), float(m.group(2))
            rows.append((e_mev * 1000.0, h_pSv_cm2))  # MeV -> keV

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("# Photon ambient dose equivalent H*(10) per fluence, ICRP Publication 74 / ICRU Report 57\n")
        f.write(f"# source: {URL}\n")
        f.write("# energy_keV,h_star_10_pSv_cm2\n")
        for e_keV, h in rows:
            f.write(f"{e_keV:g},{h:g}\n")
    print(f"{len(rows)}件を書き出しました: {OUT}")


if __name__ == "__main__":
    main()
