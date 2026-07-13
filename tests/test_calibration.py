"""線源のフォトン数校正（mAs）のテスト。"""
from __future__ import annotations

import numpy as np
import pytest

from vivemonte.scene import validate_scene
from vivemonte.transport import photon_count_through_field, run_transport

_BASE_SOURCE = {
    "kvp": 120, "position": [0, -180, 140], "direction": [0, 1, 0],
    "field": {"size_cm": [35, 43], "sid_cm": 180}, "filtration_mm_al": 2.5,
}
_BASE_GEOMETRY = [{
    "name": "slab", "shape": "box", "material": "water",
    "center": [0, 0, 140], "size_cm": [20, 20, 20],
}]


def test_photon_count_scales_linearly_with_mas():
    src1 = {**_BASE_SOURCE, "mas": 1.0}
    src10 = {**_BASE_SOURCE, "mas": 10.0}
    n1 = photon_count_through_field(src1)
    n10 = photon_count_through_field(src10)
    assert np.isclose(n10, n1 * 10.0, rtol=1e-6)


def test_photon_count_scales_with_field_area():
    src_small = {**_BASE_SOURCE, "mas": 1.0,
                 "field": {"size_cm": [10, 10], "sid_cm": 180}}
    src_large = {**_BASE_SOURCE, "mas": 1.0,
                 "field": {"size_cm": [20, 20], "sid_cm": 180}}
    n_small = photon_count_through_field(src_small)
    n_large = photon_count_through_field(src_large)
    assert np.isclose(n_large, n_small * 4.0, rtol=1e-6)


def test_photon_count_requires_mas():
    src = {**_BASE_SOURCE}
    with pytest.raises(ValueError):
        photon_count_through_field(src)


def test_scene_validates_mas_field():
    raw = {"source": {**_BASE_SOURCE, "mas": 5.0}, "geometry": _BASE_GEOMETRY}
    scene = validate_scene(raw)
    assert scene.ok
    assert scene.raw["source"]["mas"] == 5.0


def test_scene_rejects_negative_mas():
    raw = {"source": {**_BASE_SOURCE, "mas": -1.0}, "geometry": _BASE_GEOMETRY}
    scene = validate_scene(raw)
    assert not scene.ok


def test_run_transport_without_mas_has_no_calibration():
    raw = {"source": _BASE_SOURCE, "geometry": _BASE_GEOMETRY}
    scene = validate_scene(raw)
    result = run_transport(scene, n_histories=5_000, seed=1)
    assert result.n_photons_real is None


def test_run_transport_with_mas_calibrates():
    raw = {"source": {**_BASE_SOURCE, "mas": 4.0}, "geometry": _BASE_GEOMETRY}
    scene = validate_scene(raw)
    result = run_transport(scene, n_histories=5_000, seed=1, dose_grid=True, grid_resolution_cm=10.0)
    assert result.n_photons_real is not None
    assert result.n_photons_real > 0
    assert np.isclose(result.n_photons_real, photon_count_through_field(scene.raw["source"]))

    scale = result.n_photons_real / result.n_histories
    total_kerma_per_history_MeV = result.grid.total_kerma_MeV() / result.n_histories
    total_kerma_absolute_MeV = total_kerma_per_history_MeV * scale
    assert total_kerma_absolute_MeV > 0
    assert np.isfinite(total_kerma_absolute_MeV)
