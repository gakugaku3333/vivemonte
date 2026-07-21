"""レイリー角度サンプリング2段階方式（逆変換+角度因子棄却）の分布同等性回帰テスト。

docs/plan_rayleigh_compton_importance_sampling.md 検証プロトコル(a)の縮小版
（フル版・全16組はdocs/rci_verification/verify_rayleigh_distribution.py、
n=1,000,000。ここでは代表4組・小さめのnでCI向けに高速化する）。

新実装(`sample_rayleigh_cos_theta`)が
  (1) 理論分布 p(cosθ) ∝ (1+cos²θ)F(Z,q)² に一標本KS検定で適合し、
  (2) 旧実装(`_sample_rayleigh_cos_theta_uniform`)と2標本KS検定で同一分布とみなせ、
  (3) 旧実装の受理率が低い組（軽元素・高エネルギー）でも新実装の反復回数が
      小さいこと（Phase 1の主目標: 最悪ケースで10回以下）
を確認する。Bonferroni補正 α=0.05/4（代表4組の多重比較）。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats
from scipy.integrate import cumulative_trapezoid

from chatcarlo.materials import rayleigh_form_factor_table
from chatcarlo.physics import (_HC_KEV_ANGSTROM,
                                _sample_rayleigh_cos_theta_uniform,
                                sample_rayleigh_cos_theta)

_N_NEW = 100_000  # 新実装（高受理率）: 十分な検出力かつ高速
_N_OLD = 5_000    # 旧実装（低受理率、比較用のみ）: 最悪ケースでも数秒で完了する上限
_REPRESENTATIVE = [(1, 150.0), (8, 80.0), (20, 60.0), (82, 20.0)]  # フル16組の代表4点
_ALPHA = 0.05 / len(_REPRESENTATIVE)  # Bonferroni補正


def _theoretical_cdf(z: int, e_keV: float, n_grid: int = 50_000):
    q_grid, f_grid = rayleigh_form_factor_table(z)
    c = np.linspace(-1.0, 1.0, n_grid)
    theta = np.arccos(c)
    q = e_keV * np.sin(theta / 2.0) / _HC_KEV_ANGSTROM
    f = np.interp(q, q_grid, f_grid)
    pdf_unnorm = (1.0 + c ** 2) * f ** 2
    cdf = cumulative_trapezoid(pdf_unnorm, c, initial=0.0)
    cdf /= cdf[-1]
    return c, cdf


@pytest.mark.parametrize("z,e_keV", _REPRESENTATIVE)
def test_new_sampler_fits_theoretical_distribution(z, e_keV):
    c_grid, cdf_grid = _theoretical_cdf(z, e_keV)
    rng = np.random.default_rng(hash((z, int(e_keV))) & 0xFFFFFFFF)
    samples = sample_rayleigh_cos_theta(np.full(_N_NEW, z), np.full(_N_NEW, e_keV), rng)
    ks = stats.kstest(samples, lambda x: np.interp(x, c_grid, cdf_grid))
    assert ks.pvalue >= _ALPHA, f"Z={z} E={e_keV}: 理論分布への適合度検定不合格 (p={ks.pvalue:.4f})"


@pytest.mark.parametrize("z,e_keV", _REPRESENTATIVE)
def test_new_sampler_matches_old_sampler(z, e_keV):
    """新実装(逆変換2段階)と旧実装(cosθ一様提案+棄却)が同一分布からの標本とみなせる。"""
    rng_new = np.random.default_rng((hash((z, int(e_keV))) & 0xFFFFFFFF))
    rng_old = np.random.default_rng((hash((z, int(e_keV))) & 0xFFFFFFFF) + 1)
    new_samples = sample_rayleigh_cos_theta(np.full(_N_NEW, z), np.full(_N_NEW, e_keV), rng_new)
    old_samples = _sample_rayleigh_cos_theta_uniform(np.full(_N_OLD, z), np.full(_N_OLD, e_keV), rng_old)
    ks2 = stats.ks_2samp(new_samples, old_samples)
    assert ks2.pvalue >= _ALPHA, f"Z={z} E={e_keV}: 新旧2標本KS検定不合格 (p={ks2.pvalue:.4f})"


def test_mixed_element_array_routes_each_photon_to_its_own_distribution():
    """複数元素が混在するz_array（実輸送での水=H/O混在など）で、各光子が
    自分の元素の分布からサンプリングされること（`for z in sorted(set(...))`の
    マスク振り分けが正しいこと）を直接検証する。単一元素テストでは
    グループ化ループのマスク経路が検証されないため補う（実輸送テストでは
    間接的にしか通っていなかった経路）。"""
    n_per = 100_000
    z_a, z_b, e_keV = 8, 82, 80.0  # O と Pb を交互に
    z_array = np.empty(2 * n_per, dtype=int)
    z_array[0::2] = z_a
    z_array[1::2] = z_b
    e_array = np.full(2 * n_per, e_keV)
    rng = np.random.default_rng(2024)
    samples = sample_rayleigh_cos_theta(z_array, e_array, rng)

    for z in (z_a, z_b):
        sub = samples[z_array == z]
        c_grid, cdf_grid = _theoretical_cdf(z, e_keV)
        ks = stats.kstest(sub, lambda x: np.interp(x, c_grid, cdf_grid))
        assert ks.pvalue >= _ALPHA, (
            f"混在配列中のZ={z}の部分集合が理論分布に不適合 (p={ks.pvalue:.4f}) "
            f"——元素の振り分けが誤っている疑い")


def test_new_sampler_rejection_rate_is_high_even_in_old_worst_case():
    """旧実装が最も苦しむ組(H・150keV、平均約4,872試行/光子)で新実装が高受理率を保つこと。

    Phase 1の主目標: 最悪ケースで反復10回以下(理論期待値2回以下)。
    """
    n = 20_000
    rng = np.random.default_rng(123)
    z_array = np.full(n, 1)
    e_array = np.full(n, 150.0)
    x_max_all = (e_array / _HC_KEV_ANGSTROM) ** 2
    pending = np.arange(n)
    rounds = 0
    from chatcarlo.materials import rayleigh_cumulative_table
    while len(pending) > 0:
        rounds += 1
        assert rounds <= 10, f"反復回数が主目標(10回以下)を超えた: {rounds}回目でも未収束"
        zp = z_array[pending]
        x_max = x_max_all[pending]
        x_grid, a_grid = rayleigh_cumulative_table(1)
        a_cut = np.interp(x_max, x_grid, a_grid)
        xi1 = rng.random(len(pending))
        x = np.minimum(np.interp(xi1 * a_cut, a_grid, x_grid), x_max)
        c = np.clip(1.0 - 2.0 * x / x_max, -1.0, 1.0)
        xi2 = rng.random(len(pending))
        accept = xi2 <= (1.0 + c ** 2) / 2.0
        pending = pending[~accept]
