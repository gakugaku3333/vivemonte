"""ChatCarlo側の相互検証実行スクリプト（150 keV鉛筆ビーム＋水スラブ10cm）。

幾何・線質・許容基準はPREREGISTRATION.md参照。Phase 1（60 keV水10cm）と同一
パターン・同一厚さ、エネルギーのみ変更。実行:

    PYTHONPATH=. .venv/bin/python docs/egs5_crosscheck/water150kev/run_chatcarlo_water150.py
"""
from __future__ import annotations

import math

import numpy as np

from chatcarlo.geometry import Geometry
from chatcarlo.materials import linear_mu
from chatcarlo.transport import transport_photons

MATERIAL, THICKNESS_CM, ENERGY_KEV, N, SEED = "water", 10.0, 150.0, 500_000, 1


def main() -> None:
    # bbox_margin_cmを既定の50cmのまま使うと、スラブ脱出後も世界境界まで
    # 追加で~50cmの空気中を飛行することになり、その空気減弱がEGS5の
    # tutor5パターン（スラブ外は真空）とずれる（20keV版で約5.5%減弱と判明、
    # docs/lessons_learned.md参照）。真空境界に近づけるため、bbox_marginと
    # 線源ギャップを意図的に極小にする。
    geom = Geometry([{
        "name": "slab", "shape": "box", "material": MATERIAL,
        "center": [0.0, 0.0, 0.0],
        "size_cm": [THICKNESS_CM, 100.0, 100.0],
    }], bbox_margin_cm=0.01)
    rng = np.random.default_rng(SEED)
    pos = np.tile(np.array([-THICKNESS_CM / 2 - 0.01, 0.0, 0.0]), (N, 1))
    dirv = np.tile(np.array([1.0, 0.0, 0.0]), (N, 1))
    energy = np.full(N, ENERGY_KEV)

    result = transport_photons(pos, dirv, energy, geom, rng)
    uncollided = np.sum(result.escaped & (result.n_scatter == 0)) / N

    mu = float(linear_mu(MATERIAL, np.array([ENERGY_KEV]))[0])
    expected = math.exp(-mu * THICKNESS_CM)
    stderr = math.sqrt(uncollided * (1 - uncollided) / N)

    print(f"mu(xraylib, water, {ENERGY_KEV} keV) = {mu:.6f} /cm")
    print(f"解析解 Beer-Lambert exp(-mu*t)      = {expected * 100:.4f}%")
    print(f"ChatCarlo MC 一次透過率 (n={N}, seed={SEED}) = {uncollided * 100:.4f}%"
          f"  (二項近似の統計誤差 ±{stderr * 100:.4f}pp)")


if __name__ == "__main__":
    main()
