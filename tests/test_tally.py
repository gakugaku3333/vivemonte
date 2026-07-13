"""ボクセル線量タリーのテスト。

1. accumulate_track_length単体: 区間全体のエネルギーが（どのボクセルに
   割り振られようと）過不足なく積算されることを厳密に検証。
2. 統合テスト: 同一の輸送シミュレーション中で、相互作用点ごとに集計する
   衝突推定量（transport_photons.energy_deposited、既存の物理検証済み）と、
   飛程積分によるカーマtrack-length estimator（グリッド積算）が、
   同じ物理量（全カーマ）を異なる手法で推定していることを利用し、
   統計誤差の範囲で一致することを確認する（独立した交差検証）。
"""
from __future__ import annotations

import numpy as np

from vivemonte.geometry import Geometry
from vivemonte.tally import VoxelGrid, accumulate_track_length
from vivemonte.transport import transport_photons


def test_accumulate_exact_energy_conservation():
    grid = VoxelGrid.from_bbox(np.array([0.0, 0.0, 0.0]), np.array([10.0, 10.0, 10.0]), resolution_cm=1.0)
    n = 500
    origin = np.tile(np.array([0.5, 0.5, 0.0]), (n, 1))
    direction = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    length_cm = np.full(n, 7.3)          # ボクセル境界をまたぐ長さ
    weight = np.full(n, 12.0)            # keV/cm

    accumulate_track_length(grid, origin, direction, length_cm, weight, substep_cm=0.3)

    expected_total = np.sum(length_cm * weight)
    assert np.isclose(grid.kerma_keV.sum(), expected_total, rtol=1e-9)


def test_accumulate_partial_out_of_grid():
    grid = VoxelGrid.from_bbox(np.array([0.0, 0.0, 0.0]), np.array([10.0, 10.0, 10.0]), resolution_cm=1.0)
    n = 200
    # グリッド外(x=-5)からグリッドを突き抜けて外(x=15)へ抜ける長い区間
    origin = np.tile(np.array([-5.0, 5.0, 5.0]), (n, 1))
    direction = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    length_cm = np.full(n, 20.0)
    weight = np.full(n, 5.0)

    accumulate_track_length(grid, origin, direction, length_cm, weight, substep_cm=0.5)

    total_in_grid = np.sum(length_cm * weight) * (10.0 / 20.0)  # グリッド内は全長の半分
    assert grid.kerma_keV.sum() > 0
    assert grid.kerma_keV.sum() < np.sum(length_cm * weight)
    assert np.isclose(grid.kerma_keV.sum(), total_in_grid, rtol=0.05)


def test_grid_kerma_matches_collision_estimator():
    """同一の輸送で、track-length estimator(グリッド)とcollision estimator
    (相互作用点ごとの直接集計)が同じ全カーマを別手法で推定し、統計誤差内で一致する。"""
    material, thickness, energy_keV, n = "water", 15.0, 60.0, 300_000
    geom = Geometry([{
        "name": "slab", "shape": "box", "material": material,
        "center": [0.0, 0.0, 0.0], "size_cm": [thickness, 60.0, 60.0],
    }])
    grid = VoxelGrid.from_bbox(geom.bbox_min, geom.bbox_max, resolution_cm=3.0)

    rng = np.random.default_rng(7)
    pos = np.tile(np.array([-thickness / 2 - 5.0, 0.0, 0.0]), (n, 1))
    dirv = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    energy = np.full(n, energy_keV)

    result = transport_photons(pos, dirv, energy, geom, rng, grid=grid)

    collision_total_keV = sum(result.energy_deposited.values())
    tracklength_total_keV = grid.kerma_keV.sum()

    rel_diff = abs(tracklength_total_keV - collision_total_keV) / collision_total_keV
    assert rel_diff < 0.05
