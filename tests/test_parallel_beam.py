"""parallel照射野（平行ビーム近似、非発散）の検証。"""
from __future__ import annotations

import numpy as np
import pytest

from vivemonte.scene import field_corners, validate_scene
from vivemonte.source import sample_source_photons

_PAR_SRC = {
    "spectrum": [{"energy_keV": 60.0, "weight": 1.0}],
    "position": [0.0, 0.0, 15.0],
    "direction": [0.0, 1.0, 0.0],
    "field": {"shape": "parallel", "size_cm": [10.0, 10.0]},
}


def test_parallel_directions_are_all_identical():
    """非発散ビーム: 全光子の方向が完全に一致する（cone/rectのような広がりがない）。"""
    rng = np.random.default_rng(1)
    pos, dirs, energy = sample_source_photons(_PAR_SRC, 5000, rng)
    assert np.allclose(dirs, [0.0, 1.0, 0.0])
    assert np.all(energy == 60.0)


def test_parallel_origins_uniform_over_size_cm():
    """光子の始点は position を中心にsize_cmの面上に一様分布し、ビーム軸(y)は不変。"""
    rng = np.random.default_rng(2)
    n = 100_000
    pos, dirs, energy = sample_source_photons(_PAR_SRC, n, rng)
    assert np.allclose(pos[:, 1], 0.0)  # 全光子がy=0(position)面から出発
    assert pos[:, 0].min() >= -5.0 - 1e-9 and pos[:, 0].max() <= 5.0 + 1e-9
    assert pos[:, 2].min() >= 10.0 - 1e-9 and pos[:, 2].max() <= 20.0 + 1e-9
    se = (10.0 / np.sqrt(12.0)) / np.sqrt(n)
    assert abs(pos[:, 0].mean() - 0.0) < 5 * se
    assert abs(pos[:, 2].mean() - 15.0) < 5 * se


def test_field_corners_returns_position_plane_for_parallel():
    pts = field_corners(_PAR_SRC)
    assert len(pts) == 4
    for p in pts:
        assert p[1] == pytest.approx(0.0)  # SIDを持たず、position自体の面
        assert abs(p[0]) == pytest.approx(5.0)
        assert abs(p[2] - 15.0) == pytest.approx(5.0)


def test_scene_validation_accepts_parallel_without_sid():
    scene = validate_scene({
        "source": {**_PAR_SRC},
        "geometry": [{"name": "p", "shape": "box", "material": "water",
                      "size_cm": [30, 20, 30], "center": [0, 10, 15]}],
    })
    assert scene.ok, scene.errors
