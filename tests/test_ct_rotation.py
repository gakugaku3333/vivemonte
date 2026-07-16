"""CTガントリー回転の離散近似（source.rotation）の検証。

sample_source_photonsが、isocenter周りの円周上に焦点位置を離散的に
配置し、ビーム方向が常にisocenter側を向くことを確認する。
"""
from __future__ import annotations

import numpy as np

from chatcarlo.source import sample_source_photons

_ISO = np.array([0.0, 0.0, 110.0])
_SRC = {
    "kvp": 120.0,
    "filtration_mm_al": 2.5,
    "position": [0.0, 0.0, 170.0],
    "direction": [0.0, 0.0, -1.0],
    "rotation": {"isocenter": _ISO.tolist(), "axis": "y", "n_angles": 36},
    "field": {"size_cm": [50.0, 4.0], "sid_cm": 60.0},
}


def test_rotating_source_positions_lie_on_isocenter_circle():
    rng = np.random.default_rng(7)
    pos, dirv, energy = sample_source_photons(_SRC, 5000, rng)

    radius_expected = 60.0  # |position - isocenter|
    # 回転軸(y)方向にはisocenterからずれない
    assert np.allclose(pos[:, 1], _ISO[1], atol=1e-9)
    r = np.linalg.norm(pos[:, [0, 2]] - _ISO[[0, 2]], axis=1)
    assert np.allclose(r, radius_expected, atol=1e-6)


def test_rotating_source_covers_discrete_angles():
    rng = np.random.default_rng(11)
    n_angles = 36
    pos, dirv, energy = sample_source_photons({**_SRC, "rotation": {**_SRC["rotation"], "n_angles": n_angles}},
                                               20000, rng)
    theta = np.arctan2(pos[:, 0] - _ISO[0], pos[:, 2] - _ISO[2])
    bin_idx = np.round(theta / (2 * np.pi / n_angles)).astype(int) % n_angles
    n_distinct = len(set(bin_idx.tolist()))
    # 十分な光子数があれば36門すべてが埋まるはず
    assert n_distinct == n_angles


def test_rotating_source_direction_points_toward_isocenter_axis():
    rng = np.random.default_rng(13)
    pos, dirv, energy = sample_source_photons(_SRC, 2000, rng)
    # フィールド中心（su=sv=0相当）ではないので厳密一致は求めないが、
    # 焦点からisocenter方向への正射影が支配的（前方散乱寄りの発散ビーム）であるはず
    to_iso = _ISO[None, :] - pos
    to_iso /= np.linalg.norm(to_iso, axis=1, keepdims=True)
    cos_angle = np.sum(dirv * to_iso, axis=1)
    assert np.all(cos_angle > 0.9)


def test_no_rotation_key_falls_back_to_static_source():
    rng = np.random.default_rng(17)
    src = {k: v for k, v in _SRC.items() if k != "rotation"}
    pos, dirv, energy = sample_source_photons(src, 100, rng)
    assert np.allclose(pos, np.asarray(src["position"]))


def test_continuous_rotation_when_n_angles_omitted():
    """n_angles省略時は連続一様: 細かいビンで刻んでもすべて埋まり、離散クラスタがない。"""
    rng = np.random.default_rng(19)
    src = {**_SRC, "rotation": {"isocenter": _ISO.tolist(), "axis": "y"}}
    pos, dirv, energy = sample_source_photons(src, 50000, rng)
    theta = np.arctan2(pos[:, 0] - _ISO[0], pos[:, 2] - _ISO[2])
    # 離散36門なら720ビン中36ビンしか埋まらない。連続なら（統計的に）全ビン埋まる
    hist, _ = np.histogram(theta, bins=720, range=(-np.pi, np.pi))
    assert np.all(hist > 0)
    # 半径は変わらず円周上
    r = np.linalg.norm(pos[:, [0, 2]] - _ISO[[0, 2]], axis=1)
    assert np.allclose(r, 60.0, atol=1e-6)


def test_helical_scan_translates_uniformly_along_axis():
    """scan_length_cm指定時（ヘリカル位相平均近似）: 回転軸方向に一様分布で平行移動。"""
    rng = np.random.default_rng(23)
    scan = 30.0
    src = {**_SRC, "rotation": {"isocenter": _ISO.tolist(), "axis": "y",
                                 "scan_length_cm": scan}}
    n = 50000
    pos, dirv, energy = sample_source_photons(src, n, rng)
    y = pos[:, 1] - _ISO[1]
    assert y.min() >= -scan / 2 and y.max() <= scan / 2
    # 一様分布の平均0・標準偏差 scan/sqrt(12) に一致（5σ許容）
    se_mean = (scan / np.sqrt(12.0)) / np.sqrt(n)
    assert abs(y.mean()) < 5 * se_mean
    assert abs(y.std() - scan / np.sqrt(12.0)) < 0.2
    # 回転面内の半径はスキャン移動の影響を受けない
    r = np.linalg.norm(pos[:, [0, 2]] - _ISO[[0, 2]], axis=1)
    assert np.allclose(r, 60.0, atol=1e-6)
