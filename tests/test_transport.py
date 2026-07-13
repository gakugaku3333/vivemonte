"""輸送カーネルの解析解ベンチマーク — Beer-Lambert則との照合。

単色鉛筆ビームが既知厚みのスラブを無相互作用で通過する割合（一次透過率）は
exp(-mu*t) に一致するはずである。これは transport_photons の τ 消費ロジックが
正しいかどうかの直接的な検証になる。
"""
from __future__ import annotations

import numpy as np

from vivemonte.geometry import Geometry
from vivemonte.materials import linear_mu
from vivemonte.transport import transport_photons


def _pencil_beam_slab(material: str, thickness_cm: float, energy_keV: float, n: int, seed: int):
    geom = Geometry([{
        "name": "slab", "shape": "box", "material": material,
        "center": [0.0, 0.0, 0.0],
        "size_cm": [thickness_cm, 100.0, 100.0],
    }])
    rng = np.random.default_rng(seed)
    pos = np.tile(np.array([-thickness_cm / 2 - 10.0, 0.0, 0.0]), (n, 1))
    dirv = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    energy = np.full(n, energy_keV)
    result = transport_photons(pos, dirv, energy, geom, rng)
    return result


def test_beer_lambert_water_60kev():
    material, thickness, energy_keV, n = "water", 10.0, 60.0, 200_000
    result = _pencil_beam_slab(material, thickness, energy_keV, n, seed=1)
    uncollided = np.sum(result.escaped & (result.n_scatter == 0)) / n

    mu = linear_mu(material, energy_keV)
    expected = np.exp(-mu * thickness)

    stderr = np.sqrt(expected * (1 - expected) / n)
    assert abs(uncollided - expected) < 5 * stderr


def test_beer_lambert_lead_80kev():
    material, thickness, energy_keV, n = "lead", 0.2, 80.0, 200_000
    result = _pencil_beam_slab(material, thickness, energy_keV, n, seed=2)
    uncollided = np.sum(result.escaped & (result.n_scatter == 0)) / n

    mu = linear_mu(material, energy_keV)
    expected = np.exp(-mu * thickness)

    stderr = np.sqrt(expected * (1 - expected) / n)
    assert abs(uncollided - expected) < 5 * stderr


def test_energy_conservation():
    """吸収エネルギー総量は入射エネルギー総量を超えない（散逸のみ、生成なし）。"""
    material, thickness, energy_keV, n = "bone", 5.0, 50.0, 50_000
    result = _pencil_beam_slab(material, thickness, energy_keV, n, seed=3)
    total_in = n * energy_keV
    total_deposited = sum(result.energy_deposited.values())
    assert 0 <= total_deposited <= total_in
