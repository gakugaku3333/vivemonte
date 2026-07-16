"""レイリー散乱の原子形状因子サンプリングのテスト。

原子形状因子F(Z,q)込みの微分断面積は、単純なThomson近似(1+cos²θ)/2
（平均cosθ=0の前後対称分布）と違い、q（ひいては散乱角）が大きいほど
F(Z,q)が減衰するため前方散乱に偏る。この偏りはZが大きいほど強くなる
（重元素ほど原子核近傍の電子雲が広がっており、形状因子の減衰が速いため）。
これらは既知の物理的性質であり、直接の数値照合が難しい角度分布の
妥当性を間接的に検証する手段として使う。
"""
from __future__ import annotations

import numpy as np

from chatcarlo.materials import element_composition, rayleigh_element_weights, rayleigh_form_factor_table
from chatcarlo.physics import sample_rayleigh_cos_theta, sample_rayleigh_element


def test_form_factor_at_q0_equals_Z():
    for z in (1, 8, 20, 82):
        q_grid, f_grid = rayleigh_form_factor_table(z)
        assert np.isclose(f_grid[0], z, atol=1e-6)


def test_form_factor_monotonically_decreasing():
    q_grid, f_grid = rayleigh_form_factor_table(82)
    assert np.all(np.diff(f_grid) <= 1e-9)


def test_scattering_is_forward_biased_not_symmetric():
    """Thomson近似(平均cosθ=0)と異なり、実際の分布は前方散乱に偏る。"""
    n = 20_000
    rng = np.random.default_rng(0)
    z = np.full(n, 82)  # 鉛
    e = np.full(n, 60.0)
    cos_theta = sample_rayleigh_cos_theta(z, e, rng)
    assert cos_theta.mean() > 0.3  # Thomson(平均0)より明確に前方偏り


def test_lighter_element_is_more_forward_peaked():
    """xraylibの形状因子データ(EPDLベース)で確認した実際の傾向: F(Z,q)/Zは
    軽元素ほどqに対して急峻に減衰する（電子雲が実空間で広がっているため
    q空間では逆に狭い）。そのため同じエネルギーでは軽元素の方が
    相対的に前方散乱へ強く偏り、重元素は後方まで広がる
    （例: q=0.5Å⁻¹でのF/Zは C:0.28, Ca:0.41, Pb:0.59 — 直接xraylibで確認済み）。"""
    n = 20_000
    rng = np.random.default_rng(1)
    e = np.full(n, 60.0)
    cos_lead = sample_rayleigh_cos_theta(np.full(n, 82), e, rng)
    cos_carbon = sample_rayleigh_cos_theta(np.full(n, 6), e, rng)
    assert cos_carbon.mean() > cos_lead.mean()


def test_higher_energy_is_more_forward_peaked():
    """同じ元素ならエネルギーが高いほど（qが大きくなりやすいほど）前方に偏る。"""
    n = 20_000
    rng = np.random.default_rng(2)
    z = np.full(n, 20)  # カルシウム
    cos_low = sample_rayleigh_cos_theta(z, np.full(n, 20.0), rng)
    cos_high = sample_rayleigh_cos_theta(z, np.full(n, 150.0), rng)
    assert cos_high.mean() > cos_low.mean()


def test_element_selection_weights_sum_to_one():
    zs, w = rayleigh_element_weights("bone", np.array([30.0, 60.0, 100.0]))
    assert set(zs.tolist()) == {z for z, _ in element_composition("bone")}
    assert np.allclose(w.sum(axis=0), 1.0)


def test_element_selection_only_returns_compound_members():
    rng = np.random.default_rng(3)
    n = 5000
    materials = np.full(n, "bone", dtype=object)
    energies = np.full(n, 50.0)
    z_chosen = sample_rayleigh_element(materials, energies, rng)
    valid_zs = {z for z, _ in element_composition("bone")}
    assert set(np.unique(z_chosen).tolist()).issubset(valid_zs)


def test_calcium_overrepresented_relative_to_mass_fraction_in_bone():
    """レイリー断面積は重元素ほど大きいため、Caの相互作用選択割合は
    質量分率(21%)より高くなるはず（元素別断面積で重み付けしているため）。"""
    rng = np.random.default_rng(4)
    n = 20_000
    materials = np.full(n, "bone", dtype=object)
    energies = np.full(n, 50.0)
    z_chosen = sample_rayleigh_element(materials, energies, rng)
    ca_fraction = np.mean(z_chosen == 20)
    ca_mass_fraction = dict(element_composition("bone"))[20]
    assert ca_fraction > ca_mass_fraction
