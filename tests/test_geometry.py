"""Geometry.nearest_object_distance_cm のテスト。"""
from __future__ import annotations

import numpy as np

from vivemonte.geometry import Geometry


def test_nearest_object_distance_outside_box():
    geom = Geometry([{
        "name": "box", "shape": "box", "material": "water",
        "center": [0.0, 0.0, 0.0], "size_cm": [10.0, 10.0, 10.0],
    }])
    # 箱の表面(x=5)から10cm離れた点 -> 最短距離は10cm
    assert np.isclose(geom.nearest_object_distance_cm([15.0, 0.0, 0.0]), 10.0)


def test_nearest_object_distance_zero_when_inside():
    geom = Geometry([{
        "name": "box", "shape": "box", "material": "water",
        "center": [0.0, 0.0, 0.0], "size_cm": [10.0, 10.0, 10.0],
    }])
    assert geom.nearest_object_distance_cm([0.0, 0.0, 0.0]) == 0.0


def test_nearest_object_distance_picks_closest_of_multiple():
    geom = Geometry([
        {"name": "near", "shape": "sphere", "material": "water",
         "center": [10.0, 0.0, 0.0], "radius_cm": 2.0},
        {"name": "far", "shape": "sphere", "material": "water",
         "center": [-100.0, 0.0, 0.0], "radius_cm": 2.0},
    ])
    # near球の表面までの距離(10-2=8cm)がfar球より近い
    assert np.isclose(geom.nearest_object_distance_cm([0.0, 0.0, 0.0]), 8.0)


def test_nearest_object_distance_none_when_no_objects():
    geom = Geometry.__new__(Geometry)
    geom.geoms = []
    geom.background = "air"
    assert geom.nearest_object_distance_cm([0.0, 0.0, 0.0]) is None
