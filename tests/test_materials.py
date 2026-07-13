"""断面積データの検証 — NIST公表値とのスポット照合。

参照値は physics.nist.gov XAAMDI テーブル（Hubbell & Seltzer）の生値。
μ/ρ は xraylib、μen/ρ は同梱NISTテーブル経由で、両者が一次ソースと
一致することを保証する。
"""
import numpy as np
import pytest

from vivemonte.materials import linear_mu, mu_en_rho, mu_rho

# (材料, keV, NIST μ/ρ, NIST μen/ρ) — テーブルのグリッド点なので厳密一致を要求
NIST_REFERENCE = [
    ("water", 60.0, 2.059e-1, 3.190e-2),
    ("water", 100.0, 1.707e-1, 2.546e-2),
    ("aluminum", 60.0, 2.778e-1, 1.099e-1),
    ("aluminum", 100.0, 1.704e-1, 3.794e-2),
    ("lead", 100.0, 5.549e0, 1.976e0),
    ("soft_tissue", 60.0, 2.048e-1, 3.264e-2),
    ("bone", 60.0, 3.148e-1, 1.400e-1),
    ("air", 60.0, 1.875e-1, 3.041e-2),
]


# 生体組織は組成規格が供給源で異なる（xraylib=ICRP、NIST XAAMDI=ICRU-44）ため
# μ/ρ が最大約2%ずれる。これは実在する物理的差異なので許容幅を分ける。
_LOOSE = {"soft_tissue", "bone", "lung"}


@pytest.mark.parametrize("mat,e,ref_mu,ref_muen", NIST_REFERENCE)
def test_mu_rho_matches_nist(mat, e, ref_mu, ref_muen):
    rel = 0.02 if mat in _LOOSE else 0.01
    assert mu_rho(mat, e)[0] == pytest.approx(ref_mu, rel=rel)


@pytest.mark.parametrize("mat,e,ref_mu,ref_muen", NIST_REFERENCE)
def test_mu_en_rho_matches_nist(mat, e, ref_mu, ref_muen):
    # μen/ρ: 同梱テーブルのグリッド点そのものなので0.1%以内
    assert mu_en_rho(mat, e)[0] == pytest.approx(ref_muen, rel=1e-3)


def test_loglog_interpolation_between_grid_points():
    # グリッド間(70keV)の補間値が両隣の値の間に入ること
    v50, v70, v80 = (mu_en_rho("water", e)[0] for e in (50.0, 70.0, 80.0))
    assert v80 < v70 < v50


def test_water_hvl_at_60kev_sanity():
    # 60keV単色の水のHVL ≈ ln2/μ ≈ 3.37cm（教科書値）
    hvl = np.log(2) / linear_mu("water", 60.0)[0]
    assert 3.2 < hvl < 3.5


def test_unknown_material_raises_helpful_error():
    with pytest.raises(ValueError, match="候補"):
        mu_rho("unobtanium", 60.0)
