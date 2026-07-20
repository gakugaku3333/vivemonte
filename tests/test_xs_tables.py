"""断面積テーブル化（docs/plan_transport_speedup.md Phase 1）の同等性検証。

輸送カーネルは相互作用サンプリングのたびにxraylibを光子ごとに呼んでいたのを、
元素単位でエネルギー格子上に事前テーブル化しPCHIP補間で引く方式に変更した
（chatcarlo/materials.py の `_element_xs_tables` 系）。数値ロジックを
一切変えない高速化のはずなので、テーブル補間値がxraylib直接呼び出しと
（吸収端の極小近傍を除き）十分一致することをここで担保する。

実装過程で判明した落とし穴（このテストで再発防止する）:
1. xraylib.EdgeEnergy()の公称値と、実際にCS_Photoがジャンプする位置は完全には
   一致しない（Al K端で最大約0.6eV[相対0.038%]、Zn L端で約0.2%のズレを実測）。
2. 亜鉛のL3端(1.02keV)とL2端(1.04keV)の間では、光電断面積がジャンプ直後に
   単純減衰へは向かわず約1.8倍まで滑らかに上昇する「近接閾値構造」があり、
   基本格子密度だけでは補間しきれず最大約2%の誤差を生んだ。
どちらも `_locate_photo_jump_kev`（粗スキャン+二分探索）と`_NEAR_THRESHOLD_*`
（ジャンプ直後の追加補強点）で解消済み。
"""
import numpy as np
import pytest
import xraylib

from chatcarlo.materials import (MATERIAL_ALIASES, _XAAMDI_FILES,
                                  _element_cs, _element_edge_energies_kev,
                                  _XS_GRID_MAX_KEV, _XS_GRID_MIN_KEV,
                                  element_composition, mu_rho, resolve)

_EDGE_EXCLUDE_REL = 5e-3  # 吸収端±0.5%は「意図的な曖昧窓」として除外（現物窓は±0.001%）
_TOLERANCE_REL = 2e-3     # 実測worst caseは0.06〜0.07%だが、将来材料追加への余裕を見て0.2%とする


def _registered_element_zs() -> list[int]:
    """MATERIAL_ALIASES/_XAAMDI_FILESに登録済みの全材料が使う元素番号の集合。"""
    zs = set()
    for name in sorted(set(MATERIAL_ALIASES) | set(_XAAMDI_FILES)):
        try:
            zs |= {z for z, _ in element_composition(name)}
        except ValueError:
            continue  # resolve不能な登録名（既存の別問題、本テストの対象外）
    return sorted(zs)


_ELEMENT_ZS = _registered_element_zs()


@pytest.fixture(scope="module")
def dense_energy_sample():
    rng = np.random.default_rng(20260720)
    return np.exp(rng.uniform(np.log(_XS_GRID_MIN_KEV), np.log(_XS_GRID_MAX_KEV), 5000))


@pytest.mark.parametrize("z", _ELEMENT_ZS)
@pytest.mark.parametrize("kind,direct_fn", [
    ("photo", xraylib.CS_Photo), ("compt", xraylib.CS_Compt), ("rayl", xraylib.CS_Rayl),
])
def test_element_xs_table_matches_xraylib_away_from_edges(z, kind, direct_fn, dense_energy_sample):
    edges = _element_edge_energies_kev(z)

    def near_any_edge(e):
        return any(abs(e - ed) / ed < _EDGE_EXCLUDE_REL for ed in edges)

    mask = np.array([not near_any_edge(e) for e in dense_energy_sample])
    e_test = dense_energy_sample[mask]
    if len(e_test) == 0:
        pytest.skip("この元素の吸収端密度が高すぎてテスト点が残らなかった")
    direct = np.array([direct_fn(z, e) for e in e_test])
    table = _element_cs(z, e_test, kind)
    reldiff = np.abs(table - direct) / np.where(direct > 0, direct, 1.0)
    assert reldiff.max() < _TOLERANCE_REL, (
        f"Z={z} kind={kind}: 最大相対誤差{reldiff.max()*100:.4f}% "
        f"at E={e_test[np.argmax(reldiff)]:.4f} keV")


@pytest.mark.parametrize("material", sorted(set(MATERIAL_ALIASES)))
def test_compound_mu_rho_matches_xraylib_cp_away_from_edges(material, dense_energy_sample):
    """化合物のmu_rho（混合則ベース）が、xraylibの*_CP直接呼び出しと一致すること。"""
    name, _, is_elem = resolve(material)
    if is_elem:
        direct_fn = lambda e: xraylib.CS_Total(xraylib.SymbolToAtomicNumber(name), e)
    else:
        direct_fn = lambda e: xraylib.CS_Total_CP(name, e)

    all_edges = []
    for z, _ in element_composition(material):
        all_edges += list(_element_edge_energies_kev(z))

    def near_any_edge(e):
        return any(abs(e - ed) / ed < _EDGE_EXCLUDE_REL for ed in all_edges)

    e_test = dense_energy_sample[:1500]
    mask = np.array([not near_any_edge(e) for e in e_test])
    e_clean = e_test[mask]
    direct = np.array([direct_fn(e) for e in e_clean])
    table = mu_rho(material, e_clean)
    reldiff = np.abs(table - direct) / np.where(direct > 0, direct, 1.0)
    assert reldiff.max() < _TOLERANCE_REL, (
        f"{material}: 最大相対誤差{reldiff.max()*100:.4f}% at E={e_clean[np.argmax(reldiff)]:.4f} keV")


def test_lead_k_edge_jump_is_sharp_not_smeared():
    """鉛K端(88.0keV、遮蔽計算で最重要)の不連続がテーブル補間でなだらかに
    均されていないこと（Phase 1着手時に実際に起きたバグの回帰テスト:
    エネルギー幅を広く取った素朴な補間だと端を跨いで直線的に均してしまい、
    最大62%もの誤差が出た）。
    """
    z = 82
    below = xraylib.CS_Photo(z, 87.9)
    above = xraylib.CS_Photo(z, 88.1)
    table_below = _element_cs(z, np.array([87.9]), "photo")[0]
    table_above = _element_cs(z, np.array([88.1]), "photo")[0]
    assert table_below == pytest.approx(below, rel=_TOLERANCE_REL)
    assert table_above == pytest.approx(above, rel=_TOLERANCE_REL)
    assert above / below > 4.0  # 実際に大きな不連続であること自体の確認


def test_zinc_near_threshold_structure_is_captured():
    """亜鉛L3-L2端間の近接閾値構造（ジャンプ後に約1.8倍まで滑らかに上昇する
    実在の物理挙動）が、基本格子密度だけでは捉えきれず最大約2%の誤差を
    生んでいたバグの回帰テスト。"""
    z = 30
    for e in (1.025, 1.030, 1.035, 1.040):
        direct = xraylib.CS_Photo(z, e)
        table = _element_cs(z, np.array([e]), "photo")[0]
        assert table == pytest.approx(direct, rel=_TOLERANCE_REL), f"E={e}"
