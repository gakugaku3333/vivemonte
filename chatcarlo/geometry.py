"""解析ジオメトリー — レイ・プリミティブ交差とロケーション判定。

ボクセル化はしない。box/cylinder/sphere は閉形式の交差計算で扱い、
輸送カーネルは光子ごとに「次の境界までの距離」を都度求めて進む
（解析面トラッキング）。重なりは scene.yaml のリスト順で後勝ち
（後方の物体が前方を上書き）とし、リストに含まれない領域は
background（既定 air）とみなす。

座標系は scene.py と同じ: cm単位、z軸が鉛直上向き。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

_AXES = {"x": 0, "y": 1, "z": 2}


def _intersect_box(o, d, center, size, eps=1e-12):
    n = o.shape[0]
    c = np.asarray(center, dtype=float)
    s = np.asarray(size, dtype=float)
    lo = c - s / 2.0
    hi = c + s / 2.0
    t_enter = np.full(n, -np.inf)
    t_exit = np.full(n, np.inf)
    miss = np.zeros(n, dtype=bool)
    for k in range(3):
        dk = d[:, k]
        ok = o[:, k]
        parallel = np.abs(dk) < eps
        miss |= parallel & ((ok < lo[k]) | (ok > hi[k]))
        dk_safe = np.where(parallel, 1.0, dk)
        ta = (lo[k] - ok) / dk_safe
        tb = (hi[k] - ok) / dk_safe
        tmin_k = np.where(parallel, -np.inf, np.minimum(ta, tb))
        tmax_k = np.where(parallel, np.inf, np.maximum(ta, tb))
        t_enter = np.maximum(t_enter, tmin_k)
        t_exit = np.minimum(t_exit, tmax_k)
    hit = (~miss) & (t_enter <= t_exit)
    return np.where(hit, t_enter, np.nan), np.where(hit, t_exit, np.nan)


def _intersect_sphere(o, d, center, radius):
    oc = o - np.asarray(center, dtype=float)
    b = np.einsum("ij,ij->i", oc, d)
    c = np.einsum("ij,ij->i", oc, oc) - radius ** 2
    disc = b * b - c
    hit = disc >= 0
    sq = np.sqrt(np.where(hit, disc, 0.0))
    t1, t2 = -b - sq, -b + sq
    return (np.where(hit, np.minimum(t1, t2), np.nan),
            np.where(hit, np.maximum(t1, t2), np.nan))


def _intersect_cylinder(o, d, center, radius, height, axis, eps=1e-12):
    ai = _AXES[axis]
    oi = [k for k in range(3) if k != ai]
    c = np.asarray(center, dtype=float)
    n = o.shape[0]

    ox = o[:, oi[0]] - c[oi[0]]
    oy = o[:, oi[1]] - c[oi[1]]
    dx = d[:, oi[0]]
    dy = d[:, oi[1]]
    a = dx * dx + dy * dy
    b = ox * dx + oy * dy
    cc = ox * ox + oy * oy - radius ** 2

    near_zero_a = a < eps
    disc = b * b - a * cc
    hit_r = (~near_zero_a) & (disc >= 0)
    sq = np.sqrt(np.where(hit_r, disc, 0.0))
    a_safe = np.where(near_zero_a, 1.0, a)
    t1 = (-b - sq) / a_safe
    t2 = (-b + sq) / a_safe
    t_enter_r = np.where(hit_r, np.minimum(t1, t2), -np.inf)
    t_exit_r = np.where(hit_r, np.maximum(t1, t2), np.inf)

    inside_axis_line = near_zero_a & (cc <= 0)
    t_enter_r = np.where(inside_axis_line, -np.inf, t_enter_r)
    t_exit_r = np.where(inside_axis_line, np.inf, t_exit_r)
    miss_r = (~hit_r) & (~inside_axis_line)

    lo_h, hi_h = c[ai] - height / 2.0, c[ai] + height / 2.0
    oz, dz = o[:, ai], d[:, ai]
    parallel_z = np.abs(dz) < eps
    miss_z = parallel_z & ((oz < lo_h) | (oz > hi_h))
    dz_safe = np.where(parallel_z, 1.0, dz)
    ta = (lo_h - oz) / dz_safe
    tb = (hi_h - oz) / dz_safe
    t_enter_z = np.where(parallel_z, -np.inf, np.minimum(ta, tb))
    t_exit_z = np.where(parallel_z, np.inf, np.maximum(ta, tb))

    t_enter = np.maximum(t_enter_r, t_enter_z)
    t_exit = np.minimum(t_exit_r, t_exit_z)
    miss = miss_r | miss_z | (t_enter > t_exit)
    return np.where(miss, np.nan, t_enter), np.where(miss, np.nan, t_exit)


def _contains_box(points, center, size):
    c = np.asarray(center, dtype=float)
    s = np.asarray(size, dtype=float)
    lo, hi = c - s / 2, c + s / 2
    return np.all((points >= lo) & (points <= hi), axis=1)


def _contains_sphere(points, center, radius):
    return np.sum((points - np.asarray(center, dtype=float)) ** 2, axis=1) <= radius ** 2


def _contains_cylinder(points, center, radius, height, axis):
    ai = _AXES[axis]
    oi = [k for k in range(3) if k != ai]
    c = np.asarray(center, dtype=float)
    r2 = (points[:, oi[0]] - c[oi[0]]) ** 2 + (points[:, oi[1]] - c[oi[1]]) ** 2
    return (r2 <= radius ** 2) & (np.abs(points[:, ai] - c[ai]) <= height / 2)


@dataclass
class Geometry:
    """scene.yaml の geometry リストを解析ジオメトリーとして扱うラッパー。

    material_at: 点が属する材料を判定（リスト後方が優先、既定は background）
    next_boundary: 光子が次に材料境界を跨ぐまでの距離（世界境界＝脱出も含む）
    """
    geoms: list[dict]
    background: str = "air"
    bbox_margin_cm: float = 50.0
    bbox_min: np.ndarray = field(init=False)
    bbox_max: np.ndarray = field(init=False)

    def __post_init__(self):
        self.bbox_min, self.bbox_max = self._compute_bbox()

    def _shape_bbox(self, g):
        c = np.asarray(g["center"], dtype=float)
        shape = g["shape"]
        if shape == "box":
            half = np.asarray(g["size_cm"], dtype=float) / 2
            return c - half, c + half
        if shape == "sphere":
            r = g["radius_cm"]
            return c - r, c + r
        if shape == "cylinder":
            ai = _AXES[g.get("axis", "z")]
            r, h = g["radius_cm"], g["height_cm"]
            lo, hi = c.copy(), c.copy()
            for k in range(3):
                if k == ai:
                    lo[k] -= h / 2
                    hi[k] += h / 2
                else:
                    lo[k] -= r
                    hi[k] += r
            return lo, hi
        raise ValueError(f"未知のshape: {shape}")

    def _compute_bbox(self):
        los, his = zip(*(self._shape_bbox(g) for g in self.geoms))
        lo = np.min(np.stack(los), axis=0)
        hi = np.max(np.stack(his), axis=0)
        return lo - self.bbox_margin_cm, hi + self.bbox_margin_cm

    def _contains(self, g, points):
        shape = g["shape"]
        if shape == "box":
            return _contains_box(points, g["center"], g["size_cm"])
        if shape == "sphere":
            return _contains_sphere(points, g["center"], g["radius_cm"])
        if shape == "cylinder":
            return _contains_cylinder(points, g["center"], g["radius_cm"],
                                       g["height_cm"], g.get("axis", "z"))
        raise ValueError(f"未知のshape: {shape}")

    def _intersect(self, g, o, d):
        shape = g["shape"]
        if shape == "box":
            return _intersect_box(o, d, g["center"], g["size_cm"])
        if shape == "sphere":
            return _intersect_sphere(o, d, g["center"], g["radius_cm"])
        if shape == "cylinder":
            return _intersect_cylinder(o, d, g["center"], g["radius_cm"],
                                        g["height_cm"], g.get("axis", "z"))
        raise ValueError(f"未知のshape: {shape}")

    def nearest_object_distance_cm(self, point) -> float | None:
        """pointから、宣言済み各物体の境界ボックスまでの最短距離[cm]。物体が無ければNone。

        線源近傍のタリー値（点線源モデルの1/r²発散）が物理的に意味を持つ範囲かどうかの
        判定に使う: 線源からこの距離より近い領域には、シーン内のどの物体も存在しない
        （＝人や検出器が実在し得ない非物理的な近傍である）。
        """
        if not self.geoms:
            return None
        pos = np.asarray(point, dtype=float)
        dists = []
        for g in self.geoms:
            lo, hi = self._shape_bbox(g)
            closest = np.clip(pos, lo, hi)
            dists.append(float(np.linalg.norm(pos - closest)))
        return min(dists)

    def material_at(self, points: np.ndarray) -> np.ndarray:
        """点(N,3) -> 材料名の配列(N,)。リスト後方が優先、既定は background。"""
        mat = np.full(points.shape[0], self.background, dtype=object)
        for g in self.geoms:
            m = self._contains(g, points)
            if np.any(m):
                mat[m] = g["material"]
        return mat

    def next_boundary(self, o: np.ndarray, d: np.ndarray, eps: float = 1e-6):
        """光子(N)ごとに次の境界までの距離と、それが世界脱出かどうかを返す。

        戻り値: (t_boundary(N,), is_escape(N,) bool)
        """
        n = o.shape[0]
        t_obj = np.full(n, np.inf)
        for g in self.geoms:
            t1, t2 = self._intersect(g, o, d)
            for t in (t1, t2):
                valid = ~np.isnan(t) & (t > eps)
                t_obj = np.where(valid & (t < t_obj), t, t_obj)
        center = (self.bbox_min + self.bbox_max) / 2
        size = self.bbox_max - self.bbox_min
        _, t_exit_world = _intersect_box(o, d, center, size)
        t_exit_world = np.where(np.isnan(t_exit_world), np.inf, t_exit_world)
        escape = t_exit_world <= t_obj
        return np.minimum(t_obj, t_exit_world), escape
