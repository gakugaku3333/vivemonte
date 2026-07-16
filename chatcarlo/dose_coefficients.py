"""外部被ばく防護のための線量換算係数 — 周辺線量当量 H*(10)。

出典: ICRP Publication 74 (1996) / ICRU Report 57 (1998)。
テーブル値は OpenMC プロジェクト（MITライセンス）が同梱する
`openmc/data/dose/icrp74/photons_H10.txt` から取得した（値自体は
ICRP74公表値の転記。取得元・手順は scripts/fetch_h_star_10.py 参照）。

H*(10)はカーマとは異なる量（フルエンスベースの防護量）であり、
`chatcarlo/tally.py` の吸収線量タリーとは別に、フルエンスの
track-length estimatorとして積算する（chatcarlo/transport.py参照）。

診断X線領域ではH*(10)は個人線量当量Hp(10,0°)・実効線量E(AP)と
数値的にほぼ一致することが知られている（Otto, JINST 2019,
https://arxiv.org/abs/1906.05411 — Eph <= 6 MeVで成立）。
"""
from __future__ import annotations

import functools
from pathlib import Path

import numpy as np

_DATA_PATH = Path(__file__).resolve().parent / "data" / "h_star_10" / "photons_icrp74.csv"


@functools.lru_cache(maxsize=1)
def _table() -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(_DATA_PATH, delimiter=",", comments="#")
    return data[:, 0], data[:, 1]  # E_keV, h*(10)/Φ [pSv cm²]


def h_star_10_per_fluence(energies_keV) -> np.ndarray:
    """周辺線量当量換算係数 h*(10)/Φ [pSv・cm²]（log-log補間）。

    テーブル範囲は10〜10000 keV。範囲外は端点値でクランプする
    （診断領域はフィルタ後の低エネルギー端が10 keV未満になっても
    寄与が無視できるほど小さいため、この近似で十分）。
    """
    e_tab, h_tab = _table()
    e = np.atleast_1d(np.asarray(energies_keV, dtype=float))
    e_clamped = np.clip(e, e_tab[0], e_tab[-1])
    return np.exp(np.interp(np.log(e_clamped), np.log(e_tab), np.log(h_tab)))
