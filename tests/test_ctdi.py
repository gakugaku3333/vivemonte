"""CTDIvol校正（chatcarlo/ctdi.py）の検証。

物理サニティ（head>body、周辺>中心）と、mAs+SpekPyフルエンスという
独立な絶対校正経路との桁の相互突き合わせを行う。
「符号・有限性だけのテストは桁バグを検出できない」（docs/lessons_learned.md）。
"""
from __future__ import annotations

import functools

import numpy as np
import pytest

from chatcarlo.ctdi import ctdi_per_history_Gy, effective_histories_from_ctdi
from chatcarlo.scene import validate_scene

_SRC = {
    "kvp": 120.0,
    "filtration_mm_al": 2.5,
    "anode_angle_deg": 12.0,
    "position": [0.0, 0.0, 170.0],
    "direction": [0.0, 0.0, -1.0],
    "rotation": {"isocenter": [0.0, 0.0, 110.0], "axis": "y", "scan_length_cm": 30.0},
    "field": {"size_cm": [50.0, 4.0], "sid_cm": 60.0},
}


@functools.lru_cache(maxsize=None)
def _ctdi(phantom: str):
    return ctdi_per_history_Gy(_SRC, phantom=phantom, n_histories=80_000, seed=31)


def test_ctdi_positive_and_finite():
    ctdiw, d_c, d_p = _ctdi("body")
    assert np.isfinite(ctdiw) and ctdiw > 0
    assert np.isfinite(d_c) and d_c > 0
    assert np.isfinite(d_p) and d_p > 0


def test_head_phantom_higher_dose_than_body():
    """Ø16cmファントムはØ32cmより減弱が少なく、同一線源なら線量が高い。"""
    ctdiw_body, _, _ = _ctdi("body")
    ctdiw_head, _, _ = _ctdi("head")
    assert ctdiw_head > ctdiw_body


def test_periphery_exceeds_center_in_body_phantom():
    """ボウタイなし・Ø32cm・120kVでは表面下1cmの周辺孔が中心孔より高線量。"""
    _, d_center, d_periph = _ctdi("body")
    assert d_periph > d_center


def test_effective_histories_consistency():
    """実効光子数×CTDIw/history ≈ 入力CTDIvol（換算の自己整合性）。"""
    ctdi_vol_mGy = 12.0
    src = {**_SRC, "ctdi_vol_mGy": ctdi_vol_mGy, "ctdi_phantom": "body"}
    n_eff = effective_histories_from_ctdi(src, seed=31)
    ctdiw, _, _ = ctdi_per_history_Gy(src, phantom="body", n_histories=200_000, seed=31)
    assert n_eff * ctdiw * 1e3 == pytest.approx(ctdi_vol_mGy, rel=0.15)


def test_ctdi_order_of_magnitude_against_spekpy_mas_path():
    """独立な絶対校正経路（mAs+SpekPyフルエンス）とCTDI経路の桁が一致すること。

    総mAs=300・コリメーション4cm・スキャン30cm → 実効mAs=300×4/30=40。
    120kV bodyファントムのCTDIvol/実効mAsの文献代表値 ~0.05-0.2 mGy/mAs
    → 予測CTDIvol ~2-8 mGy。ボウタイなし等のモデル差を含めても
    0.5〜50 mGyの帯を外れたら桁バグを疑う。
    """
    pytest.importorskip("spekpy")
    from chatcarlo.source import photon_count_through_field

    ctdiw_per_history, _, _ = _ctdi("body")
    n_photons = photon_count_through_field({**_SRC, "mas": 300.0})
    predicted_ctdivol_mGy = ctdiw_per_history * n_photons * 1e3
    assert 0.5 < predicted_ctdivol_mGy < 50.0, (
        f"予測CTDIvol={predicted_ctdivol_mGy:.3g} mGy — mAs経路とCTDI経路の"
        "桁が不整合。スケーリングバグ（n_histories二重除算等）を疑うこと"
    )


def test_scene_validation_rejects_mas_and_ctdi_together():
    raw = {
        "source": {**_SRC, "mas": 300.0, "ctdi_vol_mGy": 12.0,
                   "position": [0, 0, 170], "direction": [0, 0, -1]},
        "geometry": [{"name": "p", "shape": "sphere", "material": "water",
                      "radius_cm": 10, "center": [0, 0, 110]}],
    }
    scene = validate_scene(raw)
    assert not scene.ok
    assert any("同時に指定できません" in str(e) for e in scene.errors)


def test_scene_validation_requires_rotation_for_ctdi():
    src = {k: v for k, v in _SRC.items() if k != "rotation"}
    raw = {
        "source": {**src, "ctdi_vol_mGy": 12.0},
        "geometry": [{"name": "p", "shape": "sphere", "material": "water",
                      "radius_cm": 10, "center": [0, 0, 110]}],
    }
    scene = validate_scene(raw)
    assert not scene.ok
    assert any("rotation" in str(e) for e in scene.errors)
