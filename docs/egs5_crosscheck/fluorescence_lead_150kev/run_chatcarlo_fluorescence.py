"""K殻蛍光X線エネルギー依存性EGS5相互検証（鉛、150 keV）: ChatCarlo側の実行スクリプト。

幾何・線質・許容基準はPREREGISTRATION.md参照。150 keV単色鉛筆ビーム->
鉛スラブ(厚さ0.135cm)、蛍光ON/OFFの2条件、n=1,000,000、seed=1。
transport_photonsを直接叩く。docs/egs5_crosscheck/fluorescence/(100 keV版)
と同一パターン（エネルギー・厚さのみ変更）。

統計誤差はdocs/lessons_learned.md「統計誤差は二項近似ではなく標本標準偏差から
直接計算する」の教訓に従い、透過エネルギー割合はhistory毎の値の標本標準偏差
(x.std(ddof=1)/sqrt(N))から直接計算する（タングステン100 keV版で後から
判明した教訓を最初から反映）。
"""
from __future__ import annotations

import json

import numpy as np

from chatcarlo.geometry import Geometry
from chatcarlo.materials import linear_mu
from chatcarlo.transport import transport_photons

ENERGY_KEV = 150.0
THICKNESS_CM = 0.135
N = 1_000_000
SEED = 1
PEAK_BAND = (72.0, 86.0)


def _slab_arrays(n, seed):
    # bbox_margin_cmを既定の50cmのまま使うと、スラブ脱出後も世界境界まで
    # 追加で~50cmの空気中を飛行することになる（water150kev検証で確認済み、
    # docs/lessons_learned.md参照）。EGS5のtutor5パターン（スラブ外は真空）
    # に合わせるため、bbox_marginと線源ギャップを意図的に極小にする。
    geometry = Geometry([{
        "name": "slab", "shape": "box", "material": "lead",
        "center": [0.0, 0.0, 0.0],
        "size_cm": [THICKNESS_CM, 100.0, 100.0],
    }], bbox_margin_cm=0.01)
    rng = np.random.default_rng(seed)
    pos = np.tile(np.array([-THICKNESS_CM / 2 - 0.001, 0.0, 0.0]), (n, 1))
    dirv = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    energy = np.full(n, ENERGY_KEV)
    return pos, dirv, energy, geometry, rng


def run(fluorescence_enabled: bool) -> dict:
    pos, dirv, energy, geometry, rng = _slab_arrays(N, SEED)
    result = transport_photons(pos, dirv, energy, geometry, rng,
                                fluorescence_enabled=fluorescence_enabled)
    e_escaped_per_history = np.where(result.escaped, result.final_energy, 0.0)
    e_escaped = float(np.sum(e_escaped_per_history))
    e_in_total = N * ENERGY_KEV

    escaped_e = result.final_energy[result.escaped]
    in_peak = np.mean((escaped_e >= PEAK_BAND[0]) & (escaped_e <= PEAK_BAND[1]))
    n_escaped = int(np.sum(result.escaped))
    n_uncollided = int(np.sum(result.escaped & (result.n_scatter == 0)))

    transmitted_energy_fraction = e_escaped / e_in_total
    x_fraction = e_escaped_per_history / ENERGY_KEV
    stderr_transmitted_energy_fraction = float(x_fraction.std(ddof=1) / np.sqrt(N))

    return {
        "fluorescence_enabled": fluorescence_enabled,
        "n_histories": N,
        "seed": SEED,
        "energy_keV": ENERGY_KEV,
        "thickness_cm": THICKNESS_CM,
        "linear_mu_per_cm": float(linear_mu("lead", ENERGY_KEV)[0]),
        "n_fluorescence": int(result.n_fluorescence),
        "n_escaped": n_escaped,
        "uncollided_transmission": n_uncollided / N,
        "uncollided_transmission_stderr": float(np.sqrt(
            (n_uncollided / N) * (1 - n_uncollided / N) / N)),
        "transmitted_energy_fraction": transmitted_energy_fraction,
        "transmitted_energy_fraction_stderr": stderr_transmitted_energy_fraction,
        "peak_band_fraction_of_escaped": float(in_peak),
    }


if __name__ == "__main__":
    results = {}
    for fluor in (False, True):
        label = "on" if fluor else "off"
        r = run(fluor)
        results[label] = r
        print(f"--- fluorescence={label} ---")
        for k, v in r.items():
            print(f"  {k}: {v}")

    delta = (results["on"]["transmitted_energy_fraction"]
              - results["off"]["transmitted_energy_fraction"])
    results["delta_transmitted_energy_fraction_on_minus_off"] = delta
    print(f"\nΔ(透過エネルギー割合, ON-OFF) = {delta:.6f}")

    with open("docs/egs5_crosscheck/fluorescence_lead_150kev/chatcarlo_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print("\n書き出し: docs/egs5_crosscheck/fluorescence_lead_150kev/chatcarlo_results.json")
