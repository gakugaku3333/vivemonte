"""cone照射野（立体角一様のコーンビーム）の検証。"""
from __future__ import annotations

import numpy as np
import pytest

from vivemonte.scene import field_corners, validate_scene
from vivemonte.source import cone_half_angle_rad, sample_source_photons

_CONE_SRC = {
    "kvp": 100.0,
    "filtration_mm_al": 2.5,
    "position": [0.0, 0.0, 170.0],
    "direction": [0.0, 0.0, -1.0],
    "field": {"shape": "cone", "diameter_cm": 40.0, "sid_cm": 60.0},
}


def test_cone_directions_within_half_angle_and_unit_norm():
    rng = np.random.default_rng(41)
    pos, dirs, energy = sample_source_photons(_CONE_SRC, 20000, rng)
    assert np.allclose(np.linalg.norm(dirs, axis=1), 1.0, atol=1e-12)
    cos_half = np.cos(cone_half_angle_rad(_CONE_SRC["field"]))
    cos_theta = dirs @ np.array([0.0, 0.0, -1.0])
    assert np.all(cos_theta >= cos_half - 1e-12)


def test_cone_directions_uniform_in_solid_angle():
    """立体角一様 ⇔ cosθが[cos半頂角, 1]で一様。平均と分散を照合する。"""
    rng = np.random.default_rng(43)
    n = 100_000
    pos, dirs, energy = sample_source_photons(_CONE_SRC, n, rng)
    cos_half = np.cos(cone_half_angle_rad(_CONE_SRC["field"]))
    cos_theta = dirs @ np.array([0.0, 0.0, -1.0])
    width = 1.0 - cos_half
    se = (width / np.sqrt(12.0)) / np.sqrt(n)
    assert abs(cos_theta.mean() - (1.0 + cos_half) / 2.0) < 5 * se
    assert abs(cos_theta.std() - width / np.sqrt(12.0)) < 5 * se
    # 方位角の一様性: 横方向成分の平均はゼロ
    assert abs(dirs[:, 0].mean()) < 5e-3
    assert abs(dirs[:, 1].mean()) < 5e-3


def test_cone_composes_with_rotation():
    """cone照射野＋ガントリー回転: 各光子の方向が回転後のビーム軸（焦点→isocenter）
    から半頂角以内に収まる。"""
    rng = np.random.default_rng(47)
    iso = np.array([0.0, 0.0, 110.0])
    src = {**_CONE_SRC,
           "rotation": {"isocenter": iso.tolist(), "axis": "y"}}
    pos, dirs, energy = sample_source_photons(src, 20000, rng)
    beam_axis = iso[None, :] - pos
    beam_axis /= np.linalg.norm(beam_axis, axis=1, keepdims=True)
    cos_theta = np.sum(dirs * beam_axis, axis=1)
    cos_half = np.cos(cone_half_angle_rad(src["field"]))
    assert np.all(cos_theta >= cos_half - 1e-9)


def test_field_corners_returns_ring_for_cone():
    src = {**_CONE_SRC}
    pts = field_corners(src)
    assert len(pts) == 16
    # 全点が開口円周上（中心軸からの距離=半径20cm、SID面z=110）
    for p in pts:
        assert p[2] == pytest.approx(110.0)
        assert np.hypot(p[0], p[1]) == pytest.approx(20.0)


def _minimal_scene(field: dict) -> dict:
    return {
        "source": {"type": "xray_tube", "kvp": 100, "position": [0, 0, 170],
                   "direction": [0, 0, -1], "field": field},
        "geometry": [{"name": "p", "shape": "sphere", "material": "water",
                      "radius_cm": 10, "center": [0, 0, 110]}],
    }


def test_scene_validation_cone_requires_diameter():
    scene = validate_scene(_minimal_scene({"shape": "cone", "sid_cm": 60}))
    assert not scene.ok
    assert any("diameter_cm" in str(e) for e in scene.errors)


def test_scene_validation_rejects_unknown_field_shape():
    scene = validate_scene(_minimal_scene({"shape": "fan", "sid_cm": 60,
                                           "size_cm": [40, 40]}))
    assert not scene.ok
    assert any("field.shape" in str(e) for e in scene.errors)


def test_scene_validation_accepts_cone():
    scene = validate_scene(_minimal_scene({"shape": "cone", "diameter_cm": 40,
                                           "sid_cm": 60}))
    assert scene.ok


def test_photon_count_equal_area_cone_matches_rect():
    """mAs校正の照射野面積: 同面積のcone/rectで実光子数が一致する。"""
    pytest.importorskip("spekpy")
    from vivemonte.source import photon_count_through_field

    base = {"kvp": 100.0, "filtration_mm_al": 2.5, "mas": 10.0}
    # 面積を一致させる: π·r² = w·h → r=10cm, w=h=√(100π)
    side = float(np.sqrt(100.0 * np.pi))
    n_cone = photon_count_through_field({**base,
        "field": {"shape": "cone", "diameter_cm": 20.0, "sid_cm": 100.0}})
    n_rect = photon_count_through_field({**base,
        "field": {"shape": "rect", "size_cm": [side, side], "sid_cm": 100.0}})
    assert n_cone == pytest.approx(n_rect, rel=1e-9)
