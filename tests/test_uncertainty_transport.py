"""統計不確かさの輸送組み込み（transport.py, Phase 2）のテスト。

docs/plan_statistical_uncertainty.md Phase 2 の受入基準:
  1. 統計ON/OFFでtotal（kerma_keV/h10_track_pSv_cm3/energy_deposited_MeV等）が
     ビット一致すること（最重要・絶対制約）。
  2. Rが1/√Nでスケールすること（統計的、複数回反復平均で検証）。
  3. バッチ分割を変えてもRが同オーダーであること（統計的）。
  4. workers=1とworkers>=2でRが同オーダーであること。
  5. batch_size=1（各バッチ=1history）でのRが、実際の輸送物理から得た
     per-history寄与列のブルートフォースstd(ddof=1)と一致すること
     （Phase 1で検証済みの推定器が、実物理を通しても正しく配線されている確認）。
  6. dose_grid=Falseでも材料別付与エネルギーのSEM/Rが得られること（副目標）。

チェックポイント個別のRは統計量なのでこのファイルの多くのテストは緩めの許容誤差
（実測でキャリブレーション済み）を使う——ここは物理検証ではなく配線の正しさの
確認であることに注意。
"""
from __future__ import annotations

import numpy as np

from chatcarlo.scene import validate_scene
from chatcarlo.tally import VoxelGrid
from chatcarlo.transport import run_transport, transport_photons

_RAW_SCENE = {
    "source": {"kvp": 100, "position": [0, -50, 0], "direction": [0, 1, 0],
               "field": {"size_cm": [30, 30], "sid_cm": 100}},
    "geometry": [
        {"name": "target", "shape": "box", "material": "water",
         "center": [0, 0, 0], "size_cm": [20, 20, 20]},
    ],
}
_SCENE = validate_scene(_RAW_SCENE)


def test_totals_bit_identical_regardless_of_track_uncertainty_serial():
    """統計機構ON/OFFで、total値（絶対制約）がビット一致する。dose_grid・
    材料別付与エネルギー・吸収/脱出割合・蛍光カウントすべてを見る。"""
    r_on = run_transport(_SCENE, n_histories=20000, seed=5, dose_grid=True,
                          grid_resolution_cm=5.0, track_uncertainty=True)
    r_off = run_transport(_SCENE, n_histories=20000, seed=5, dose_grid=True,
                           grid_resolution_cm=5.0, track_uncertainty=False)

    assert np.array_equal(r_on.grid.kerma_keV, r_off.grid.kerma_keV)
    assert np.array_equal(r_on.grid.h10_track_pSv_cm3, r_off.grid.h10_track_pSv_cm3)
    assert r_on.energy_deposited_MeV == r_off.energy_deposited_MeV
    assert r_on.fraction_absorbed == r_off.fraction_absorbed
    assert r_on.fraction_escaped == r_off.fraction_escaped
    assert r_on.mean_scatter_events == r_off.mean_scatter_events
    assert r_on.n_fluorescence == r_off.n_fluorescence
    assert r_on.n_photons_real == r_off.n_photons_real


def test_totals_bit_identical_regardless_of_track_uncertainty_parallel():
    """並列パス（workers>=2）でも同様にビット一致する
    （ワーカー側のend_batchもtotalへの書き込みを一切行わないため）。"""
    r_on = run_transport(_SCENE, n_histories=16000, seed=9, n_workers=2, dose_grid=True,
                          grid_resolution_cm=5.0, track_uncertainty=True)
    r_off = run_transport(_SCENE, n_histories=16000, seed=9, n_workers=2, dose_grid=True,
                           grid_resolution_cm=5.0, track_uncertainty=False)

    assert np.array_equal(r_on.grid.kerma_keV, r_off.grid.kerma_keV)
    assert np.array_equal(r_on.grid.h10_track_pSv_cm3, r_off.grid.h10_track_pSv_cm3)
    assert r_on.energy_deposited_MeV == r_off.energy_deposited_MeV
    assert r_on.n_fluorescence == r_off.n_fluorescence


def test_default_track_uncertainty_matches_explicit_true():
    """run_transportの既定（track_uncertainty省略）はTrue（設計判断6「標準搭載」）。"""
    r_default = run_transport(_SCENE, n_histories=5000, seed=3, dose_grid=True, grid_resolution_cm=5.0)
    r_explicit = run_transport(_SCENE, n_histories=5000, seed=3, dose_grid=True, grid_resolution_cm=5.0,
                                track_uncertainty=True)
    assert np.array_equal(r_default.grid.kerma_keV, r_explicit.grid.kerma_keV)
    assert r_default.n_batches == r_explicit.n_batches
    assert r_default.energy_deposited_rel_err.keys() == r_explicit.energy_deposited_rel_err.keys()


def test_same_seed_workers_reproducible_including_uncertainty_stats():
    """同一(seed, workers)なら、totalだけでなく統計量（sum2・n_batches）も
    2回の実行でビット一致する（並列集約の加算順序がワーカー番号順で固定されて
    いることの確認、plan_phase3_parallel.mdの既存契約の統計版）。"""
    r1 = run_transport(_SCENE, n_histories=8000, seed=123, n_workers=2, dose_grid=True,
                        grid_resolution_cm=5.0, track_uncertainty=True)
    r2 = run_transport(_SCENE, n_histories=8000, seed=123, n_workers=2, dose_grid=True,
                        grid_resolution_cm=5.0, track_uncertainty=True)
    assert np.array_equal(r1.grid.kerma_sum2, r2.grid.kerma_sum2)
    assert np.array_equal(r1.grid.h10_sum2, r2.grid.h10_sum2)
    assert np.array_equal(r1.grid.n_batches_hit, r2.grid.n_batches_hit)
    assert r1.n_batches == r2.n_batches
    assert r1.energy_deposited_rel_err == r2.energy_deposited_rel_err


def _max_voxel_index():
    """代表ボクセル（大きなnで安定させた最大カーマ位置）を固定して使う。
    毎回の小さなrunごとに独自のargmaxを取ると、位置選択自体のばらつきが
    1/√Nスケーリングの検証ノイズに上乗せされるため。"""
    r_ref = run_transport(_SCENE, n_histories=200_000, seed=999, dose_grid=True,
                           grid_resolution_cm=5.0, batch_size=2000, track_uncertainty=True)
    return np.unravel_index(int(np.argmax(r_ref.grid.kerma_keV)), r_ref.grid.kerma_keV.shape)


def test_relative_error_scales_with_inverse_sqrt_n():
    """4倍のhistoryでRがおよそ半分になる（1/√Nスケーリング）。

    単発のrunでは統計揺らぎが大きいため（実測: 単発ペアでratio 0.98〜3.4と
    広く散る）、複数repの平均で検証する——推定器自体の妥当性はPhase 1で
    厳密検証済みなので、ここでは「輸送に組み込んでも1/√N則が壊れていないか」
    という配線確認が目的。"""
    voxel = _max_voxel_index()
    n_reps = 6

    def avg_R(n, seed0):
        rs = []
        for i in range(n_reps):
            r = run_transport(_SCENE, n_histories=n, seed=seed0 + i, dose_grid=True,
                               grid_resolution_cm=5.0, batch_size=2000, track_uncertainty=True)
            rs.append(r.grid.kerma_relative_error()[voxel])
        return np.mean(rs)

    r_small = avg_R(20_000, seed0=300)
    r_large = avg_R(80_000, seed0=400)

    ratio = r_small / r_large
    assert 1.3 < ratio < 3.0, f"1/√Nスケーリングから外れている: ratio={ratio:.3f} (期待値付近2.0)"


def test_relative_error_order_of_magnitude_invariant_to_batch_size():
    """同一n_historiesでbatch_sizeだけ変えても（M=20 vs M=10相当）、
    Rの値は同オーダーに留まる（厳密一致ではなく、点推定としての妥当な範囲内）。"""
    voxel = _max_voxel_index()
    r_fine = run_transport(_SCENE, n_histories=40_000, seed=11, dose_grid=True,
                            grid_resolution_cm=5.0, batch_size=2000, track_uncertainty=True)
    r_coarse = run_transport(_SCENE, n_histories=40_000, seed=12, dose_grid=True,
                              grid_resolution_cm=5.0, batch_size=4000, track_uncertainty=True)

    assert r_fine.n_batches == 20
    assert r_coarse.n_batches == 10

    r_fine_val = r_fine.grid.kerma_relative_error()[voxel]
    r_coarse_val = r_coarse.grid.kerma_relative_error()[voxel]
    ratio = r_fine_val / r_coarse_val
    assert 0.2 < ratio < 5.0, (
        f"batch_sizeを変えるとRのオーダーが崩れている: fine(M=20)={r_fine_val:.4f}, "
        f"coarse(M=10)={r_coarse_val:.4f}")


def test_parallel_relative_error_same_order_as_serial():
    """workers=1とworkers>=2で、同一n_historiesにおけるRが同オーダーであること
    （異なる乱数ストリーム・端数バッチの丸め方が違うため厳密一致は要求しない、
    plan_phase3_parallel.mdの「統計的同等」契約の延長）。"""
    voxel = _max_voxel_index()
    r_serial = run_transport(_SCENE, n_histories=40_000, seed=7, dose_grid=True,
                              grid_resolution_cm=5.0, batch_size=2000, track_uncertainty=True,
                              n_workers=1)
    r_parallel = run_transport(_SCENE, n_histories=40_000, seed=7, dose_grid=True,
                                grid_resolution_cm=5.0, batch_size=2000, track_uncertainty=True,
                                n_workers=2)

    r_s = r_serial.grid.kerma_relative_error()[voxel]
    r_p = r_parallel.grid.kerma_relative_error()[voxel]
    assert np.isfinite(r_s) and np.isfinite(r_p)
    ratio = r_s / r_p
    assert 0.2 < ratio < 5.0, f"serial R={r_s:.4f}, parallel R={r_p:.4f}"


def test_energy_deposited_sem_available_without_dose_grid():
    """副目標: dose_grid=Falseでも材料別付与エネルギーのSEM/相対誤差が得られる
    （EGS5相互検証スクリプト群が自前のΣx・Σx²計算を捨てられる状態の確認）。"""
    result = run_transport(_SCENE, n_histories=20_000, seed=21, batch_size=2000,
                            dose_grid=False, track_uncertainty=True)
    assert result.grid is None
    assert result.n_batches == 10
    assert result.energy_deposited_rel_err
    for name, rel in result.energy_deposited_rel_err.items():
        assert np.isfinite(rel)
        assert rel >= 0.0
        assert name in result.energy_deposited_sem_MeV
        assert result.energy_deposited_sem_MeV[name] >= 0.0


def test_energy_deposited_sem_empty_when_uncertainty_disabled():
    result = run_transport(_SCENE, n_histories=5000, seed=1, dose_grid=False,
                            track_uncertainty=False)
    assert result.energy_deposited_sem_MeV == {}
    assert result.energy_deposited_rel_err == {}
    assert result.n_batches == 0


def test_relative_error_nan_when_batch_count_insufficient():
    """batch_size >= n_historiesならM=1になり、Rはnan（バッチ数不足）。"""
    result = run_transport(_SCENE, n_histories=5000, seed=2, batch_size=200_000,
                            dose_grid=False, track_uncertainty=True)
    assert result.n_batches == 1
    assert all(np.isnan(v) for v in result.energy_deposited_rel_err.values())


def test_bruteforce_matches_direct_transport_at_batch_size_one():
    """batch_size=1相当（各バッチ=1history）のときの相対誤差が、実際の輸送物理を
    1historyずつ実行して得たper-history寄与列のブルートフォース
    std(ddof=1)/√N/meanと一致することを、run_transportではなく低レベルAPI
    （transport_photons + VoxelGrid直接操作）で検証する
    （Phase 1で検証済みの推定器の代数が、実物理を通しても正しく配線されている
    ことの確認。run_transportは1historyごとに乱数ストリームを消費し、光子数
    サンプリングのオーバーヘッドがあるため低レベルAPIの方が高速）。"""
    material, thickness, energy_keV = "water", 10.0, 60.0
    from chatcarlo.geometry import Geometry
    geom = Geometry([{
        "name": "slab", "shape": "box", "material": material,
        "center": [0.0, 0.0, 0.0], "size_cm": [thickness, 60.0, 60.0],
    }])
    grid = VoxelGrid.from_bbox(geom.bbox_min, geom.bbox_max, resolution_cm=2.0, track_uncertainty=True)

    n_histories = 400
    rng = np.random.default_rng(17)
    prev_kerma = grid.kerma_keV.copy()
    deltas = []
    for _ in range(n_histories):
        pos = np.array([[-thickness / 2 - 5.0, 0.0, 0.0]])
        dirv = np.array([[1.0, 0.0, 0.0]])
        energy = np.array([energy_keV])
        transport_photons(pos, dirv, energy, geom, rng, grid=grid)
        grid.end_batch(1)
        delta = grid.kerma_keV - prev_kerma
        deltas.append(delta.copy())
        prev_kerma = grid.kerma_keV.copy()

    deltas = np.stack(deltas)  # (n_histories, nx, ny, nz)
    voxel = np.unravel_index(int(np.argmax(grid.kerma_keV)), grid.kerma_keV.shape)
    x = deltas[:, voxel[0], voxel[1], voxel[2]]
    assert x.sum() > 0  # 選んだボクセルに実際に寄与があること

    mean = x.mean()
    sem_brute = x.std(ddof=1) / np.sqrt(n_histories)
    expected_rel = sem_brute / mean

    got_rel = grid.kerma_relative_error()[voxel]
    assert grid.n_batches == n_histories
    assert np.isclose(got_rel, expected_rel, rtol=1e-9)
