"""ChatCarlo側の相互検証実行スクリプト（60 keV鉛筆ビーム＋水スラブ10cm）。

EGS5（run_water60 / run_water60_bound）と同一条件でChatCarloの一次透過率を計算する。
RESULTS.mdに記載した数値の再現用。実行:

    PYTHONPATH=. .venv/bin/python docs/egs5_crosscheck/run_chatcarlo_water60.py
"""
from __future__ import annotations

import math

import numpy as np

from chatcarlo.geometry import Geometry
from chatcarlo.materials import linear_mu
from chatcarlo.transport import transport_photons

MATERIAL, THICKNESS_CM, ENERGY_KEV, N, SEED = "water", 10.0, 60.0, 500_000, 1


def main() -> None:
    geom = Geometry([{
        "name": "slab", "shape": "box", "material": MATERIAL,
        "center": [0.0, 0.0, 0.0],
        "size_cm": [THICKNESS_CM, 100.0, 100.0],
    }])
    rng = np.random.default_rng(SEED)
    pos = np.tile(np.array([-THICKNESS_CM / 2 - 10.0, 0.0, 0.0]), (N, 1))
    dirv = np.tile(np.array([1.0, 0.0, 0.0]), (N, 1))
    energy = np.full(N, ENERGY_KEV)

    result = transport_photons(pos, dirv, energy, geom, rng)
    uncollided = np.sum(result.escaped & (result.n_scatter == 0)) / N

    mu = float(linear_mu(MATERIAL, np.array([ENERGY_KEV]))[0])
    expected = math.exp(-mu * THICKNESS_CM)
    # 二項分布近似による統計誤差（EGS5側の誤差も同じ式の手計算。EGS5の出力には
    # 統計誤差は印字されない）
    stderr = math.sqrt(uncollided * (1 - uncollided) / N)

    print(f"mu(xraylib, water, {ENERGY_KEV} keV) = {mu:.6f} /cm")
    print(f"解析解 Beer-Lambert exp(-mu*t)      = {expected * 100:.3f}%")
    print(f"ChatCarlo MC 一次透過率 (n={N}, seed={SEED}) = {uncollided * 100:.3f}%"
          f"  (二項近似の統計誤差 ±{stderr * 100:.4f}pp)")


if __name__ == "__main__":
    main()
