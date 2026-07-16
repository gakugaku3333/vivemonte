"""線源のフォトン数校正（mAs）のテスト。"""
from __future__ import annotations

import numpy as np
import pytest

from chatcarlo.scene import validate_scene
from chatcarlo.source import photon_count_through_field
from chatcarlo.transport import run_transport

_BASE_SOURCE = {
    "kvp": 120, "position": [0, -180, 140], "direction": [0, 1, 0],
    "field": {"size_cm": [35, 43], "sid_cm": 180}, "filtration_mm_al": 2.5,
}
_BASE_GEOMETRY = [{
    "name": "slab", "shape": "box", "material": "water",
    "center": [0, 0, 140], "size_cm": [20, 20, 20],
}]


def test_photon_count_scales_linearly_with_mas():
    src1 = {**_BASE_SOURCE, "mas": 1.0}
    src10 = {**_BASE_SOURCE, "mas": 10.0}
    n1 = photon_count_through_field(src1)
    n10 = photon_count_through_field(src10)
    assert np.isclose(n10, n1 * 10.0, rtol=1e-6)


def test_photon_count_scales_with_field_area():
    src_small = {**_BASE_SOURCE, "mas": 1.0,
                 "field": {"size_cm": [10, 10], "sid_cm": 180}}
    src_large = {**_BASE_SOURCE, "mas": 1.0,
                 "field": {"size_cm": [20, 20], "sid_cm": 180}}
    n_small = photon_count_through_field(src_small)
    n_large = photon_count_through_field(src_large)
    assert np.isclose(n_large, n_small * 4.0, rtol=1e-6)


def test_photon_count_requires_mas():
    src = {**_BASE_SOURCE}
    with pytest.raises(ValueError):
        photon_count_through_field(src)


def test_scene_validates_mas_field():
    raw = {"source": {**_BASE_SOURCE, "mas": 5.0}, "geometry": _BASE_GEOMETRY}
    scene = validate_scene(raw)
    assert scene.ok
    assert scene.raw["source"]["mas"] == 5.0


def test_scene_rejects_negative_mas():
    raw = {"source": {**_BASE_SOURCE, "mas": -1.0}, "geometry": _BASE_GEOMETRY}
    scene = validate_scene(raw)
    assert not scene.ok


def test_run_transport_without_mas_has_no_calibration():
    raw = {"source": _BASE_SOURCE, "geometry": _BASE_GEOMETRY}
    scene = validate_scene(raw)
    result = run_transport(scene, n_histories=5_000, seed=1)
    assert result.n_photons_real is None


def test_run_transport_with_mas_calibrates():
    raw = {"source": {**_BASE_SOURCE, "mas": 4.0}, "geometry": _BASE_GEOMETRY}
    scene = validate_scene(raw)
    result = run_transport(scene, n_histories=5_000, seed=1, dose_grid=True, grid_resolution_cm=10.0)
    assert result.n_photons_real is not None
    assert result.n_photons_real > 0
    assert np.isclose(result.n_photons_real, photon_count_through_field(scene.raw["source"]))

    # dose_per_history (=グリッド値/n_histories、1光子あたり) に n_photons_real を
    # 掛けるだけで絶対値になる。n_historiesでさらに割ってはいけない
    # （過去にここで二重に割るバグがあり、絶対値が桁で小さくなっていた。
    # 本テストは符号・有限性しか見ておらず検出できなかったため、
    # 下の桁チェックを追加した — [[lessons_learned]]参照）。
    total_kerma_per_history_MeV = result.grid.total_kerma_MeV() / result.n_histories
    total_kerma_absolute_MeV = total_kerma_per_history_MeV * result.n_photons_real
    assert total_kerma_absolute_MeV > 0
    assert np.isfinite(total_kerma_absolute_MeV)


def test_calibrated_dose_matches_spekpy_kerma_order_of_magnitude():
    """絶対値校正が正しい桁になっているかをSpekPyの自由空間カーマと突き合わせる。

    水スラブ手前（線源からの距離=170cm、スラブは減弱・散乱源になるので厳密一致は
    期待しないが、桁が合っていることを見る）でSpekPyが計算する自由空間カーマ
    （mas=4での絶対値）と、ChatCarloのボクセル線量タリー最大値を比較する。
    """
    import spekpy as sp

    from chatcarlo.geometry import Geometry
    from chatcarlo.diagnostics import dose_map_Gy

    raw = {"source": {**_BASE_SOURCE, "mas": 4.0}, "geometry": _BASE_GEOMETRY}
    scene = validate_scene(raw)
    result = run_transport(scene, n_histories=300_000, seed=1, dose_grid=True, grid_resolution_cm=1.0)

    geometry = Geometry(scene.raw["geometry"])
    dose_per_history_Gy = dose_map_Gy(result.grid, geometry) / result.n_histories
    dose_absolute_Gy = dose_per_history_Gy * result.n_photons_real

    s = sp.Spek(kvp=_BASE_SOURCE["kvp"], th=_BASE_SOURCE.get("anode_angle_deg", 12.0), z=170, mas=4.0)
    s.filter("Al", _BASE_SOURCE["filtration_mm_al"])
    expected_free_air_kerma_Gy = s.get_kerma() * 1e-6  # uGy -> Gy

    # スラブ手前面はビルドアップ・後方散乱で自由空間カーマより高くなり得るし、
    # 1cm粗さのボクセル境界には統計ノイズも乗る。ここは「桁が合っているか」
    # （factor 50以内）だけを見るゲートで、1e6ずれるようなバグを検出できれば十分。
    assert dose_absolute_Gy.max() > expected_free_air_kerma_Gy / 10
    assert dose_absolute_Gy.max() < expected_free_air_kerma_Gy * 50
