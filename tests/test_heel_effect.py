"""ヒール効果（source.heel_effect）の検証。

SpekPyの軸外計算の座標系（x<0が陽極側）は実測で確認済み:
x=-15cmでフルエンス約32%減・平均エネルギー約3.8keV硬化（120kV, th=12, z=100）。
ここでは viveMonte 側の実装が「陽極側で光子が少ない・硬い」を再現することを見る。
"""
from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("spekpy")

from vivemonte.scene import validate_scene
from vivemonte.source import sample_source_photons

# 陽極方向 = +x。SID100cm・照射野40cm角（ヒール軸±20cm）で効果がよく見える
_HEEL_SRC = {
    "kvp": 120.0,
    "filtration_mm_al": 2.5,
    "anode_angle_deg": 12.0,
    "position": [0.0, 0.0, 170.0],
    "direction": [0.0, 0.0, -1.0],
    "heel_effect": True,
    "anode_direction": [1.0, 0.0, 0.0],
    "field": {"shape": "rect", "size_cm": [40.0, 40.0], "sid_cm": 100.0},
}


def _sample(n=60_000, seed=53, **over):
    rng = np.random.default_rng(seed)
    return sample_source_photons({**_HEEL_SRC, **over}, n, rng)


def test_fewer_photons_on_anode_side():
    pos, dirs, energy = _sample()
    # 方向の+x成分が正 = SID面で陽極側(+x)に向かう光子
    n_anode = int(np.sum(dirs[:, 0] > 0))
    n_cathode = int(np.sum(dirs[:, 0] < 0))
    assert n_anode < n_cathode * 0.92, (n_anode, n_cathode)


def test_harder_spectrum_on_anode_side():
    pos, dirs, energy = _sample()
    # SID面でのヒール座標に換算して端同士を比べる（±10cmより外側）
    s = dirs[:, 0] / (-dirs[:, 2]) * 100.0
    e_anode = energy[s > 10.0].mean()
    e_cathode = energy[s < -10.0].mean()
    assert e_anode > e_cathode + 1.0, (e_anode, e_cathode)


def test_heel_off_is_symmetric():
    pos, dirs, energy = _sample(heel_effect=False)
    n_pos = int(np.sum(dirs[:, 0] > 0))
    n_neg = int(np.sum(dirs[:, 0] < 0))
    assert abs(n_pos - n_neg) < 5 * np.sqrt(len(dirs))


def test_heel_composes_with_cone_and_rotation():
    iso = [0.0, 0.0, 110.0]
    # cone Ø20 @ sid60 → 半頂角9.5度 < 陽極角12度（カットオフにかからない）
    src = {**_HEEL_SRC,
           "position": [0.0, 0.0, 170.0],
           "field": {"shape": "cone", "diameter_cm": 20.0, "sid_cm": 60.0},
           "rotation": {"isocenter": iso, "axis": "x"},  # 陽極軸(x)=回転軸と平行
           "anode_direction": [1.0, 0.0, 0.0]}
    pos, dirs, energy = _sample(n=20_000, **src)
    assert np.all(np.isfinite(dirs)) and np.all(np.isfinite(energy))
    assert np.allclose(np.linalg.norm(dirs, axis=1), 1.0, atol=1e-9)
    # 回転しても陽極側（各光子の局所ヒール軸=世界x軸: 回転軸と平行なので不変）
    # で光子が少ないまま
    n_anode = int(np.sum(dirs[:, 0] > 0))
    n_cathode = int(np.sum(dirs[:, 0] < 0))
    assert n_anode < n_cathode


def test_anode_cutoff_warns_and_still_samples():
    """照射野が陽極カットオフを超える配置では警告しつつ、抽選自体は破綻しない。"""
    src = {**_HEEL_SRC,
           "field": {"shape": "cone", "diameter_cm": 30.0, "sid_cm": 60.0}}  # 半頂角14度>12度
    with pytest.warns(UserWarning, match="陽極カットオフ"):
        pos, dirs, energy = _sample(n=5_000, **src)
    assert np.all(np.isfinite(energy))
    # カットオフ側（陽極側の外縁）にはほぼ光子が来ない
    s = dirs[:, 0] / (-dirs[:, 2] + 1e-12) * 60.0
    assert np.sum(s > 13.0) < 0.01 * len(s)


def test_scene_validation_requires_anode_direction():
    raw = {
        "source": {**{k: v for k, v in _HEEL_SRC.items() if k != "anode_direction"}},
        "geometry": [{"name": "p", "shape": "sphere", "material": "water",
                      "radius_cm": 10, "center": [0, 0, 110]}],
    }
    scene = validate_scene(raw)
    assert not scene.ok
    assert any("anode_direction" in str(e) for e in scene.errors)


def test_scene_validation_rejects_anode_parallel_to_beam():
    raw = {
        "source": {**_HEEL_SRC, "anode_direction": [0.0, 0.0, 1.0]},
        "geometry": [{"name": "p", "shape": "sphere", "material": "water",
                      "radius_cm": 10, "center": [0, 0, 110]}],
    }
    scene = validate_scene(raw)
    assert not scene.ok
    assert any("平行" in str(e) for e in scene.errors)


def test_mas_calibration_with_heel_below_central_axis_value():
    """面平均フルエンス（ヒール込み）は中心軸一様近似より小さいはず
    （中心より陽極側の低下が大きく、陰極側の増加は小さい）。"""
    from vivemonte.source import photon_count_through_field

    n_heel = photon_count_through_field({**_HEEL_SRC, "mas": 10.0})
    n_flat = photon_count_through_field(
        {**{k: v for k, v in _HEEL_SRC.items() if k != "heel_effect"}, "mas": 10.0})
    assert 0.5 * n_flat < n_heel < n_flat