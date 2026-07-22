"""バッチ統計による相対誤差推定器（tally.py）のテスト。輸送非依存（Phase 1）。

docs/plan_statistical_uncertainty.md 設計判断2の不偏推定器
    σ̂² = (Q - T²/N) / (M-1)   T=ΣS_b, Q=ΣS_b²/n_b
を検証する。3種類の検証を分けている:

1. batch_size=1（バッチ=1history）での厳密一致 — この場合のみQ=Σx_i²になり
   ブルートフォースのx.std(ddof=1)と数式的に完全一致する。
2. 不等バッチサイズでの不偏性 — batch_size=1の値とは数式的に一致しない
   （別の不偏推定量になるだけ）ため、多数回反復した平均が真の分散に収束する
   ことを統計的に確認する（厳密一致は要求しない）。
3. combine_moments — 単純加算の結合律なので厳密一致。
"""
from __future__ import annotations

import numpy as np

from chatcarlo.tally import (ScalarMoments, VoxelGrid, accumulate_track_length,
                              combine_moments, relative_error, standard_error)


def test_relative_error_matches_brute_force_at_batch_size_one():
    """batch_size=1（各バッチが1historyの寄与そのもの）ではQ=Σx_i²になり、
    relative_errorがx.std(ddof=1)/sqrt(N)/meanとrtol=1e-12で一致する。"""
    rng = np.random.default_rng(0)
    x = rng.gamma(shape=2.0, scale=3.0, size=5000)  # 非対称分布（線量寄与を模す）

    T = x.sum()
    Q = np.sum(x ** 2)  # n_b=1なのでS_b²/n_b = x_i²
    n_batches = len(x)
    n_histories = len(x)

    got = relative_error(T, Q, n_batches, n_histories)

    mean = x.mean()
    sem_brute = x.std(ddof=1) / np.sqrt(len(x))
    expected = sem_brute / mean

    assert np.isclose(got, expected, rtol=1e-12)


def test_standard_error_matches_brute_force_at_batch_size_one():
    rng = np.random.default_rng(1)
    x = rng.exponential(scale=5.0, size=3000)

    T = x.sum()
    Q = np.sum(x ** 2)
    got = standard_error(T, Q, len(x), len(x))
    expected = x.std(ddof=1) / np.sqrt(len(x))
    assert np.isclose(got, expected, rtol=1e-12)


def test_unequal_batch_sizes_unbiased_statistically():
    """不等バッチサイズのσ̂²は、batch_size=1の値とは数式的に一致しないが
    （別の不偏推定量）、多数回の反復平均は真の分散に収束するはず。"""
    true_mu, true_sigma = 10.0, 4.0
    rng = np.random.default_rng(42)
    n_reps = 4000
    var_estimates = []

    for _ in range(n_reps):
        # 不揃いなバッチサイズ: 3,1,5,2 の周期（合計11historyを1セットとして5セット分）
        batch_sizes = [3, 1, 5, 2] * 5
        n_total = sum(batch_sizes)
        x = rng.normal(true_mu, true_sigma, n_total)

        T = 0.0
        Q = 0.0
        idx = 0
        for n_b in batch_sizes:
            s_b = x[idx:idx + n_b].sum()
            T += s_b
            Q += s_b ** 2 / n_b
            idx += n_b

        var = (Q - T ** 2 / n_total) / (len(batch_sizes) - 1)
        var_estimates.append(var)

    mean_var_estimate = np.mean(var_estimates)
    assert abs(mean_var_estimate - true_sigma ** 2) / true_sigma ** 2 < 0.05


def test_combine_moments_matches_one_shot_aggregation():
    """バッチ列を2分割して個別集計→combine_momentsで合成した結果が、
    一括集計と厳密一致する（単純加算の結合律）。"""
    rng = np.random.default_rng(7)
    batch_sums = rng.gamma(2.0, 3.0, size=37)      # 37バッチ、S_b相当
    batch_ns = rng.integers(1, 6, size=37)          # 各バッチのhistory数

    def aggregate(sums, ns):
        T = float(np.sum(sums))
        Q = float(np.sum(sums ** 2 / ns))
        M = len(sums)
        N = int(np.sum(ns))
        return (T, Q, M, N)

    one_shot = aggregate(batch_sums, batch_ns)

    split = 13
    part_a = aggregate(batch_sums[:split], batch_ns[:split])
    part_b = aggregate(batch_sums[split:], batch_ns[split:])
    combined = combine_moments(part_a, part_b)

    for got, expected in zip(combined, one_shot):
        assert np.isclose(got, expected, rtol=1e-12)


def test_combine_moments_array_form():
    """T/Qがndarray（VoxelGridの場合）でもcombine_momentsが素直に加算できる。"""
    t_a = np.array([1.0, 2.0, 3.0])
    q_a = np.array([0.5, 0.5, 0.5])
    t_b = np.array([4.0, 5.0, 6.0])
    q_b = np.array([1.0, 1.0, 1.0])

    combined = combine_moments((t_a, q_a, 3, 10), (t_b, q_b, 5, 20))
    assert np.allclose(combined[0], t_a + t_b)
    assert np.allclose(combined[1], q_a + q_b)
    assert combined[2] == 8
    assert combined[3] == 30


def test_relative_error_nan_when_batches_insufficient():
    """M<2（バッチ数不足）はnan——0や無限大にせず「推定不能」と区別する。"""
    r = relative_error(T=5.0, Q=25.0, n_batches=1, n_histories=1)
    assert np.isnan(r)


def test_relative_error_nan_when_mean_is_zero():
    """寄与が真にゼロのボクセル・材料は相対誤差が定義されずnan。"""
    r = relative_error(T=0.0, Q=0.0, n_batches=5, n_histories=100)
    assert np.isnan(r)


def test_variance_clamped_nonnegative_on_rounding_noise():
    """丸め誤差でQ - T²/Nがわずかに負になっても、分散は0にクランプされ
    standard_error/relative_errorがnanや虚数にならない。"""
    n_histories = 100
    n_batches = 4
    T = 10.0
    Q = T ** 2 / n_histories - 1e-12  # わずかに負になるよう仕組む
    sem = standard_error(T, Q, n_batches, n_histories)
    assert np.isfinite(sem)
    assert sem >= 0.0


def test_relative_error_array_broadcasts_like_voxel_grid():
    """VoxelGridの使用形態（ndarray T/Q、一部ボクセルはT=0）を模した形状テスト。"""
    T = np.array([0.0, 10.0, 20.0, 0.0])
    Q = np.array([0.0, 30.0, 90.0, 0.0])
    r = relative_error(T, Q, n_batches=5, n_histories=50)
    assert r.shape == (4,)
    assert np.isnan(r[0]) and np.isnan(r[3])
    assert np.all(np.isfinite(r[[1, 2]]))


class TestScalarMoments:
    def test_accumulates_and_reports_relative_error(self):
        sm = ScalarMoments()
        rng = np.random.default_rng(3)
        per_batch = rng.gamma(2.0, 5.0, size=20)
        for s_b in per_batch:
            sm.add_batch({"water": float(s_b)}, n_histories_in_batch=1000)

        rel = sm.relative_errors()
        assert "water" in rel
        assert 0.0 < rel["water"] < 1.0

        # batch_size=1相当(各バッチのS_bが直接の観測値)なのでブルートフォースと一致
        mean = per_batch.mean()
        sem_brute = per_batch.std(ddof=1) / np.sqrt(len(per_batch))
        # ScalarMomentsのx_iは「1historyあたり」ではなく「1000historyのバッチ和」なので
        # 直接比較するにはn_batches=1history/batchの構成に揃える必要がある。ここでは
        # 別の構成(n_histories_in_batch=1)で再検証する。
        sm2 = ScalarMoments()
        for s_b in per_batch:
            sm2.add_batch({"water": float(s_b)}, n_histories_in_batch=1)
        rel2 = sm2.relative_errors()
        assert np.isclose(rel2["water"], sem_brute / mean, rtol=1e-12)

    def test_material_absent_from_some_batches_treated_as_zero(self):
        """あるバッチに現れない材料はS_b=0として扱われる（新規出現材料も同様）。"""
        sm = ScalarMoments()
        sm.add_batch({"water": 10.0, "lead": 5.0}, n_histories_in_batch=100)
        sm.add_batch({"water": 12.0}, n_histories_in_batch=100)  # leadはこのバッチで寄与ゼロ
        sm.add_batch({"water": 11.0, "lead": 6.0, "aluminum": 1.0}, n_histories_in_batch=100)

        assert sm.sums["lead"] == 5.0 + 0.0 + 6.0
        assert sm.sums["aluminum"] == 0.0 + 0.0 + 1.0
        assert sm.n_batches == 3
        assert sm.n_histories == 300

    def test_single_batch_is_nan(self):
        sm = ScalarMoments()
        sm.add_batch({"water": 42.0}, n_histories_in_batch=1000)
        rel = sm.relative_errors()
        assert np.isnan(rel["water"])

    def test_zero_histories_batch_ignored(self):
        sm = ScalarMoments()
        sm.add_batch({"water": 42.0}, n_histories_in_batch=0)
        assert sm.n_batches == 0
        assert sm.n_histories == 0


class TestVoxelGridUncertainty:
    def _grid(self, track_uncertainty):
        return VoxelGrid.from_bbox(np.array([0.0, 0.0, 0.0]), np.array([10.0, 10.0, 10.0]),
                                    resolution_cm=1.0, track_uncertainty=track_uncertainty)

    def _score_one_batch(self, grid, rng, n, weight_value):
        origin = np.tile(np.array([0.5, 0.5, 0.0]), (n, 1))
        direction = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
        length_cm = np.full(n, 3.0)
        weight = np.full(n, weight_value)
        accumulate_track_length(grid.kerma_keV, grid, origin, direction, length_cm, weight, rng)
        accumulate_track_length(grid.h10_track_pSv_cm3, grid, origin, direction, length_cm, weight, rng)

    def test_end_batch_noop_when_disabled(self):
        """track_uncertainty=False（既定）ではend_batchを呼んでも何も起きない
        （既存の呼び出し元に非侵襲、というPhase 1の非侵襲性の要件）。"""
        grid = self._grid(track_uncertainty=False)
        rng = np.random.default_rng(0)
        self._score_one_batch(grid, rng, n=100, weight_value=5.0)
        grid.end_batch(100)
        assert grid.kerma_sum2 is None
        assert grid.n_batches == 0
        assert grid.n_histories == 0

    def test_totals_bit_identical_regardless_of_tracking(self):
        """統計機構ON/OFFでtotal(kerma_keV/h10_track_pSv_cm3)がビット一致する
        （絶対制約）。end_batchはtotalへの書き込みを一切行わないスナップショット
        差分方式であることの直接確認。"""
        n_per_batch = 400
        n_batch_count = 5

        grid_off = self._grid(track_uncertainty=False)
        rng_off = np.random.default_rng(1)
        for _ in range(n_batch_count):
            self._score_one_batch(grid_off, rng_off, n_per_batch, weight_value=7.0)
            grid_off.end_batch(n_per_batch)  # track_uncertainty=Falseなのでno-op

        grid_on = self._grid(track_uncertainty=True)
        rng_on = np.random.default_rng(1)
        for _ in range(n_batch_count):
            self._score_one_batch(grid_on, rng_on, n_per_batch, weight_value=7.0)
            grid_on.end_batch(n_per_batch)

        assert np.array_equal(grid_off.kerma_keV, grid_on.kerma_keV)
        assert np.array_equal(grid_off.h10_track_pSv_cm3, grid_on.h10_track_pSv_cm3)

    def test_end_batch_normalizes_by_batch_history_count(self):
        """end_batchがn_b（バッチ内history数）で正しく正規化していることを、
        Qとは**独立な式**——分散分析の群間平方和 SSB = Σn_b(x̄_b - x̄)²——から
        検証する。

        ミューテーション試験で見つかった穴を塞ぐテスト: `end_batch`の
        `/n_histories_in_batch`を削除しても、既存テストは1件も落ちなかった。
        グリッドのRの絶対値を検証していたテストが全てn_b=1（`/n_b`が恒等変換に
        なる条件）だったため。ここではn_b>1かつバッチごとに異なるn_bを与え、
        `Q=ΣS_b²/n_b`という実装と同じ形を使わずに期待値を組む。
        """
        grid = self._grid(track_uncertainty=True)
        voxel = (3, 4, 5)
        # 各バッチの寄与S_bとhistory数n_bを決め打ちし、グリッドへ直接書き込む
        # （accumulate_track_lengthの乱数分配を挟むと期待値が解析的に書けないため、
        #  ここではend_batchの正規化だけを単離して検証する）。
        batch_sums = [12.0, 5.0, 21.0, 8.0, 15.0, 3.0]
        batch_ns = [4, 1, 7, 2, 5, 3]
        for s_b, n_b in zip(batch_sums, batch_ns):
            grid.kerma_keV[voxel] += s_b
            grid.end_batch(n_b)

        N = sum(batch_ns)
        M = len(batch_sums)
        assert grid.n_batches == M and grid.n_histories == N

        # 独立な期待値: SSB/(M-1) をバッチ平均から直接組む
        grand_mean = sum(batch_sums) / N
        ssb = sum(n_b * (s_b / n_b - grand_mean) ** 2
                  for s_b, n_b in zip(batch_sums, batch_ns))
        expected_var = ssb / (M - 1)
        expected_rel = np.sqrt(expected_var / N) / grand_mean

        assert np.isclose(grid.kerma_relative_error()[voxel], expected_rel, rtol=1e-12)

    def test_relative_error_invariant_to_batch_grouping(self):
        """同一の per-history 寄与列を、n_b=1 と n_b=K でグループ化したときの
        σ̂² が同じ母分散を推定していること（期待値として一致、n_b正規化が
        効いていなければ約√K倍ずれる）。

        `test_end_batch_normalizes_by_batch_history_count`が決め打ち値で
        代数を固定するのに対し、こちらは「グループ化しても推定対象が変わらない」
        という統計的性質の側から同じ正規化を押さえる。
        """
        rng = np.random.default_rng(5)
        n_histories, K = 600, 10
        x = rng.gamma(2.0, 3.0, size=n_histories)  # per-history寄与

        voxel = (2, 2, 2)
        grid_fine = self._grid(track_uncertainty=True)     # n_b=1
        for xi in x:
            grid_fine.kerma_keV[voxel] += xi
            grid_fine.end_batch(1)

        grid_coarse = self._grid(track_uncertainty=True)    # n_b=K
        for start in range(0, n_histories, K):
            grid_coarse.kerma_keV[voxel] += x[start:start + K].sum()
            grid_coarse.end_batch(K)

        r_fine = grid_fine.kerma_relative_error()[voxel]
        r_coarse = grid_coarse.kerma_relative_error()[voxel]
        # 同じσ²の別々の不偏推定なので厳密一致はしないが、正規化が抜けていれば
        # √K≈3.2倍ずれる。統計揺らぎ（M=60でのχ²ゆらぎ）より十分小さい帯で押さえる。
        assert 0.6 < r_fine / r_coarse < 1.7, f"fine={r_fine:.5f}, coarse={r_coarse:.5f}"

    def test_relative_error_finite_and_hit_count_tracked(self):
        grid = self._grid(track_uncertainty=True)
        rng = np.random.default_rng(2)
        for _ in range(10):
            self._score_one_batch(grid, rng, n=500, weight_value=4.0)
            grid.end_batch(500)

        assert grid.n_batches == 10
        assert grid.n_histories == 5000

        r_kerma = grid.kerma_relative_error()
        r_h10 = grid.h10_relative_error()
        hit_voxel = tuple(np.argwhere(grid.kerma_keV > 0)[0])
        assert np.isfinite(r_kerma[hit_voxel])
        assert np.isfinite(r_h10[hit_voxel])
        assert grid.n_batches_hit[hit_voxel] > 0

        untouched = tuple(np.argwhere(grid.kerma_keV == 0)[0])
        assert np.isnan(r_kerma[untouched])
        assert grid.n_batches_hit[untouched] == 0
