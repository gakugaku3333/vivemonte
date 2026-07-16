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

from chatcarlo.dose_coefficients import h_star_10_per_fluence
from chatcarlo.geometry import Geometry
from chatcarlo.tally import VoxelGrid, accumulate_track_length
from chatcarlo.transport import transport_photons


def test_accumulate_exact_energy_conservation():
    grid = VoxelGrid.from_bbox(np.array([0.0, 0.0, 0.0]), np.array([10.0, 10.0, 10.0]), resolution_cm=1.0)
    n = 500
    origin = np.tile(np.array([0.5, 0.5, 0.0]), (n, 1))
    direction = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    length_cm = np.full(n, 7.3)          # ボクセル境界をまたぐ長さ
    weight = np.full(n, 12.0)            # keV/cm

    accumulate_track_length(grid.kerma_keV, grid, origin, direction, length_cm, weight,
                             np.random.default_rng(0), substep_cm=0.3)

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

    accumulate_track_length(grid.kerma_keV, grid, origin, direction, length_cm, weight,
                             np.random.default_rng(0), substep_cm=0.5)

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


def test_h10_accumulate_exact_track_length():
    """H*(10)飛程積分も、カーマと同様にどのボクセルに割り振られようと
    区間全体の値 Σ coeff*dl を過不足なく積算する（体積正規化前の生値で検証）。"""
    grid = VoxelGrid.from_bbox(np.array([0.0, 0.0, 0.0]), np.array([10.0, 10.0, 10.0]), resolution_cm=1.0)
    n = 500
    origin = np.tile(np.array([0.5, 0.5, 0.0]), (n, 1))
    direction = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    length_cm = np.full(n, 7.3)
    coeff = h_star_10_per_fluence(60.0)[0]  # pSv・cm²
    weight = np.full(n, coeff)

    accumulate_track_length(grid.h10_track_pSv_cm3, grid, origin, direction, length_cm, weight,
                             np.random.default_rng(0), substep_cm=0.3)

    expected_total = np.sum(length_cm * weight)  # pSv・cm³
    assert np.isclose(grid.h10_track_pSv_cm3.sum(), expected_total, rtol=1e-9)


def test_h10_map_normalizes_by_voxel_volume():
    """単一ボクセル内で全長Lを直進する光子束 -> H*(10) = N * h*(10)/Φ(E) * L / V の解析解と一致。"""
    grid = VoxelGrid.from_bbox(np.array([0.0, 0.0, 0.0]), np.array([5.0, 5.0, 5.0]), resolution_cm=5.0)
    n = 1000
    length_cm = 3.0
    energy_keV = 80.0
    origin = np.tile(np.array([2.5, 2.5, 1.0]), (n, 1))
    direction = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    coeff = h_star_10_per_fluence(energy_keV)[0]
    weight = np.full(n, coeff)

    accumulate_track_length(grid.h10_track_pSv_cm3, grid, origin, direction,
                             np.full(n, length_cm), weight,
                             np.random.default_rng(0), substep_cm=0.5)

    expected_pSv = n * coeff * length_cm / grid.voxel_volume_cm3()
    assert np.isclose(grid.h10_map_pSv().sum(), expected_pSv, rtol=1e-6)


def test_boundary_start_surface_voxel_unbiased():
    """区間の始点がボクセル境界ちょうどに揃うケース（parallel照射野で全光子が
    ファントム前面から出発する状況）の回帰テスト。

    旧実装（サブステップ中点の決定的サンプリング）では量子化誤差の位相が
    全区間で同期し、表面ボクセル層が約-3%系統的に過小評価されていた
    （vive-auditor監査で発見）。層化乱数点なら不偏なので、指数分布の区間長で
    表面ボクセルへの割り当てが厳密な期待値 Σ min(L_i, 1) に統計誤差内で一致する。"""
    grid = VoxelGrid.from_bbox(np.array([0.0, 0.0, 0.0]), np.array([20.0, 10.0, 10.0]), resolution_cm=1.0)
    rng = np.random.default_rng(11)
    n = 200_000
    mfp_cm = 4.86  # 60 keV水の平均自由行程相当
    length_cm = rng.exponential(mfp_cm, n).clip(max=19.9)
    origin = np.tile(np.array([0.0, 5.5, 5.5]), (n, 1))  # 全区間がx=0境界ちょうどから出発
    direction = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    weight = np.ones(n)

    accumulate_track_length(grid.kerma_keV, grid, origin, direction, length_cm, weight,
                             np.random.default_rng(12))

    surface = grid.kerma_keV[0, 5, 5]              # x=[0,1)の表面ボクセル
    exact = np.minimum(length_cm, 1.0).sum()       # 厳密な重なり長の合計
    assert abs(surface - exact) / exact < 0.005    # 旧実装は-2.7%ずれてこのテストに落ちる


def test_run_transport_dose_grid_h10_finite_and_nonnegative():
    """実シーンでの統合テスト: H*(10)グリッドがクラッシュせず有限・非負の値を返す。"""
    from chatcarlo.scene import validate_scene
    from chatcarlo.transport import run_transport

    raw = {
        "source": {"kvp": 100, "position": [0, -50, 0], "direction": [0, 1, 0],
                    "field": {"size_cm": [30, 30], "sid_cm": 100}},
        "geometry": [
            {"name": "target", "shape": "box", "material": "water",
             "center": [0, 0, 0], "size_cm": [20, 20, 20]},
        ],
    }
    scene = validate_scene(raw)
    assert scene.ok
    result = run_transport(scene, n_histories=20_000, seed=1, dose_grid=True, grid_resolution_cm=5.0)
    h10 = result.grid.h10_map_pSv()
    assert np.all(np.isfinite(h10))
    assert np.all(h10 >= 0)
    assert h10.sum() > 0
