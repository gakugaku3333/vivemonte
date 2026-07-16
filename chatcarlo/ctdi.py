"""CTDIvol基準の絶対線量校正。

CT装置のコンソールに表示されるCTDIvol [mGy]を校正アンカーとして使う。
mAs+SpekPyフルエンス校正（photon_count_through_field）と違い、実測CTDIvolには
ボウタイフィルタ・mA変調（の平均効果）・ピッチ・実機の出力較正がすべて
折り込まれているため、CTでは汎用性・正確性ともに高い。

手順: シーンと同一の線源設定（kvp・濾過・rotation・field）で、標準CTDI
ファントム（PMMA円柱、body Ø32cm / head Ø16cm、長さ15cm、IEC 60601-2-44）を
回転中心に置いた系を内部的にシミュレートし、CTDIw相当量

    CTDIw = (1/3)·D_center + (2/3)·D_periphery

を1 historyあたり[Gy/history]で求める。D_center/D_peripheryは実測の
電離箱位置（中心軸上／表面下1cm、長さ100mm）に対応するボクセル平均線量。
実測CTDIvol[mGy]をこの値で割った比が「実効光子数」となり、mAs校正と同じ
スロット（TransportResult.n_photons_real）で per-history 値の絶対換算に使える。

近似と限界:
- 線量はdose-to-PMMA（実測は空気カーマ。μen/ρ比の差は診断領域で数%）
- ボウタイなしのため中心/周辺バランスは実機とずれる（CTDIwの重み付き平均で
  部分的に相殺されるが、患者シーン側の空間分布バイアスは残る）
- ヘリカル時のMSAD≈CTDIvol関係はスキャン範囲が測定領域＋散乱裾より
  十分長い（目安15cm以上）ことが前提
"""
from __future__ import annotations

import numpy as np

from .geometry import Geometry
from .tally import VoxelGrid

_PHANTOM_DIAMETER_CM = {"body": 32.0, "head": 16.0}
_PHANTOM_LENGTH_CM = 15.0
_CHAMBER_LENGTH_CM = 10.0   # CTDI100の100mm電離箱
_GRID_RESOLUTION_CM = 1.0


def ctdi_phantom_geometry(rot: dict, phantom: str = "body") -> Geometry:
    """回転中心に置いた標準CTDIファントム（PMMA円柱）のジオメトリー。"""
    diameter = _PHANTOM_DIAMETER_CM[phantom]
    return Geometry([{
        "name": f"ctdi_phantom_{phantom}",
        "shape": "cylinder",
        "material": "pmma",
        "axis": rot.get("axis", "z"),
        "radius_cm": diameter / 2.0,
        "height_cm": _PHANTOM_LENGTH_CM,
        "center": [float(x) for x in rot["isocenter"]],
    }])


def ctdi_per_history_Gy(src: dict, phantom: str = "body",
                         n_histories: int = 200_000,
                         seed: int | None = None) -> tuple[float, float, float]:
    """CTDIw相当量[Gy/history]と、その内訳(D_center, D_periphery)を返す。

    シーンのgeometryは使わず、線源設定srcだけを流用してファントム系を
    独立にシミュレートする（シーン内の寝台・壁等の散乱は実測CTDIにも
    含まれないので、含めないのが正しい）。
    """
    from .diagnostics import dose_map_Gy
    from .source import sample_source_photons
    from .transport import transport_photons

    rot = src.get("rotation")
    if rot is None:
        raise ValueError("CTDIvol校正には source.rotation（ガントリー回転）が必要です")

    geometry = ctdi_phantom_geometry(rot, phantom)
    radius = _PHANTOM_DIAMETER_CM[phantom] / 2.0
    iso = np.asarray(rot["isocenter"], dtype=float)
    axis_idx = {"x": 0, "y": 1, "z": 2}[rot.get("axis", "z")]
    plane = [k for k in range(3) if k != axis_idx]

    margin = 2.0
    bbox_min = iso - radius - margin
    bbox_max = iso + radius + margin
    bbox_min[axis_idx] = iso[axis_idx] - _PHANTOM_LENGTH_CM / 2.0 - margin
    bbox_max[axis_idx] = iso[axis_idx] + _PHANTOM_LENGTH_CM / 2.0 + margin
    grid = VoxelGrid.from_bbox(bbox_min, bbox_max, _GRID_RESOLUTION_CM)

    rng = np.random.default_rng(seed)
    batch = 200_000
    remaining = n_histories
    while remaining > 0:
        n = min(batch, remaining)
        remaining -= n
        pos, dirv, energy = sample_source_photons(src, n, rng)
        transport_photons(pos, dirv, energy, geometry, rng, grid=grid)

    dose = dose_map_Gy(grid, geometry) / n_histories
    centers = grid.voxel_centers().reshape(*grid.shape, 3)
    r = np.sqrt((centers[..., plane[0]] - iso[plane[0]]) ** 2
                + (centers[..., plane[1]] - iso[plane[1]]) ** 2)
    z_ok = np.abs(centers[..., axis_idx] - iso[axis_idx]) <= _CHAMBER_LENGTH_CM / 2.0

    # 中心孔: 軸から2cm以内。周辺孔: 表面下1cm（r=R-1）を挟む幅1cmの環
    center_sel = z_ok & (r < 2.0)
    periph_sel = z_ok & (r >= radius - 2.0) & (r <= radius - 1.0)
    d_center = float(dose[center_sel].mean())
    d_periph = float(dose[periph_sel].mean())
    ctdiw = d_center / 3.0 + 2.0 * d_periph / 3.0
    return ctdiw, d_center, d_periph


def effective_histories_from_ctdi(src: dict, seed: int | None = None) -> float:
    """実測CTDIvol[mGy]から実効光子数（per-history値→絶対値の換算係数）を求める。"""
    ctdi_vol_mGy = float(src["ctdi_vol_mGy"])
    phantom = src.get("ctdi_phantom", "body")
    ctdiw, _, _ = ctdi_per_history_Gy(src, phantom=phantom, seed=seed)
    return ctdi_vol_mGy * 1e-3 / ctdiw
