"""束縛コンプトン散乱（S(Z,q)非干渉性散乱関数）サンプリングのテスト。

test_rayleigh.pyと対称的な構成。S(Z,q)はq→0で0に漸近しq→∞でZに単調収束する
ため、Rayleighの形状因子F(Z,q)（q=0でZ最大、単調減少）と正反対の挙動になる:
小角散乱（小さいq）ほど強く抑制され、束縛コンプトンは自由電子KNより
後方散乱寄り・平均エネルギー移行率が大きくなる。
"""
from __future__ import annotations

import numpy as np
import pytest

from chatcarlo.materials import compton_element_weights, element_composition, incoherent_sq_table
from chatcarlo.physics import sample_compton_bound, sample_compton_element, sample_klein_nishina


def test_sq_table_is_zero_at_q0_and_approaches_Z_at_large_q():
    for z in (1, 8, 20, 82):
        q_grid, s_grid = incoherent_sq_table(z)
        assert s_grid[0] == 0.0
        assert s_grid[-1] == pytest.approx(z, rel=0.01)


def test_sq_table_monotonically_increasing():
    q_grid, s_grid = incoherent_sq_table(82)
    assert np.all(np.diff(s_grid) >= -1e-9)


def test_element_selection_weights_sum_to_one():
    zs, w = compton_element_weights("bone", np.array([30.0, 60.0, 100.0]))
    assert set(zs.tolist()) == {z for z, _ in element_composition("bone")}
    assert np.allclose(w.sum(axis=0), 1.0)


def test_element_selection_only_returns_compound_members():
    rng = np.random.default_rng(3)
    n = 5000
    materials = np.full(n, "bone", dtype=object)
    energies = np.full(n, 50.0)
    z_chosen = sample_compton_element(materials, energies, rng)
    valid_zs = {z for z, _ in element_composition("bone")}
    assert set(np.unique(z_chosen).tolist()).issubset(valid_zs)


def test_bound_compton_transfers_more_energy_than_free_electron():
    """S(q)は小角(小さいq、小さいエネルギー移行)散乱を抑制するため、
    束縛コンプトンは自由電子KNより平均エネルギー移行率<T>/Eが大きくなる
    （docs/egs5_crosscheck/check_compton_transfer.pyの机上検算と同じ結論、
    60keV水でt_KN=0.0936 < t_Sq=0.0959）。
    """
    n = 200_000
    materials = np.full(n, "water", dtype=object)
    e = np.full(n, 60.0)

    rng_bound = np.random.default_rng(0)
    eps_bound, _ = sample_compton_bound(materials, e, rng_bound)

    rng_free = np.random.default_rng(0)
    eps_free, _ = sample_klein_nishina(e, rng_free)

    t_bound = 1.0 - eps_bound.mean()
    t_free = 1.0 - eps_free.mean()
    assert t_bound > t_free


def test_bound_compton_mean_transfer_matches_independent_reference():
    """docs/egs5_crosscheck/check_compton_transfer.pyの独立机上検算
    (S(q)重み付きKN微分断面積の解析的角度積分)と統計誤差内(3sigma)で一致する
    ことを確認する。参照値は60keV水でt_Sq=0.095928（同スクリプト既報値）。
    """
    n = 2_000_000
    materials = np.full(n, "water", dtype=object)
    e = np.full(n, 60.0)
    rng = np.random.default_rng(42)
    eps, _ = sample_compton_bound(materials, e, rng)
    t = 1.0 - eps
    mean_t = t.mean()
    sem_t = t.std(ddof=1) / np.sqrt(n)
    t_ref = 0.095928
    assert abs(mean_t - t_ref) < 3 * sem_t


def test_bound_compton_more_backward_biased_than_free_electron():
    """S(q)による小角抑制の直接的な帰結として、束縛コンプトンの散乱角分布は
    自由電子KNより後方寄り（平均cosθが小さい）になる。"""
    n = 100_000
    materials = np.full(n, "water", dtype=object)
    e = np.full(n, 60.0)

    rng_bound = np.random.default_rng(7)
    _, cos_bound = sample_compton_bound(materials, e, rng_bound)

    rng_free = np.random.default_rng(7)
    _, cos_free = sample_klein_nishina(e, rng_free)

    assert cos_bound.mean() < cos_free.mean()


def test_hydrogen_overrepresented_relative_to_mass_fraction_in_bone():
    """コンプトン断面積は原子あたりの電子数(≈Z)にほぼ比例する一方、質量あたり
    では電子密度Z/Aで決まる。水素はZ/A=1で全元素中最大（重元素は中性子過剰で
    Z/A≈0.4〜0.5に下がる）ため、水素の相互作用選択割合は質量分率(4.7%)より
    明確に高くなる（Rayleigh散乱がZ^2近傍でスケールし重元素ほど過大代表される
    のとは対照的な挙動、test_rayleigh.pyのCa過大代表テストと比較）。"""
    rng = np.random.default_rng(4)
    n = 20_000
    materials = np.full(n, "bone", dtype=object)
    energies = np.full(n, 50.0)
    z_chosen = sample_compton_element(materials, energies, rng)
    h_fraction = np.mean(z_chosen == 1)
    h_mass_fraction = dict(element_composition("bone"))[1]
    assert h_fraction > h_mass_fraction
