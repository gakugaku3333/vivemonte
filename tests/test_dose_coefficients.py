"""H*(10)換算係数（ICRP74公表値）のテスト。"""
from __future__ import annotations

import numpy as np

from chatcarlo.dose_coefficients import h_star_10_per_fluence


def test_table_points_match_icrp74_exactly():
    # ICRP74 Table値そのもの（scripts/fetch_h_star_10.py で取得したCSVの検証）
    known = {10.0: 0.061, 30.0: 0.81, 60.0: 0.51, 100.0: 0.61, 150.0: 0.89}
    for e_keV, expected in known.items():
        got = h_star_10_per_fluence(e_keV)[0]
        assert np.isclose(got, expected, rtol=1e-6)


def test_has_local_minimum_near_60_kev():
    """h*(10)/Φ は10〜20keVで急峻に立ち上がり、60keV付近に極小を持つ
    （ICRU球内の光電吸収と深部到達のトレードオフによる既知の形状）。"""
    e = np.array([10.0, 20.0, 60.0, 150.0])
    h = h_star_10_per_fluence(e)
    assert h[1] > h[0]          # 10keV -> 20keVで急上昇
    assert h[2] < h[1]          # 20keV -> 60keVで低下（極小に向かう）
    assert h[3] > h[2]          # 60keV -> 150keVで再上昇


def test_out_of_range_clamped_not_error():
    low = h_star_10_per_fluence(1.0)
    high = h_star_10_per_fluence(50000.0)
    assert np.isfinite(low[0]) and np.isfinite(high[0])


def test_vectorized_matches_scalar():
    energies = np.array([10.0, 50.0, 100.0])
    batch = h_star_10_per_fluence(energies)
    singles = np.array([h_star_10_per_fluence(e)[0] for e in energies])
    assert np.allclose(batch, singles)
