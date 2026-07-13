"""ボクセル吸収線量タリー — 輸送ジオメトリーとは独立なグリッド。

輸送は解析面トラッキング（vivemonte/geometry.py, vivemonte/transport.py）
で行い、スコアリングだけをこのモジュールの一様グリッドに刻む。

track-length estimator（trackごとの飛程積分）でカーマを積算する:
  K += E * (μen/ρ) * ρ * dl
区間内はエネルギー・材料とも一定なので積分自体に離散化誤差はないが、
「その区間がどのボクセルに何cm分入っているか」の空間分配は
サブステップ分割による近似（サブステップ長を細かくするほど厳密に収束）。
電子飛程を無視するカーマ近似のため、カーマ＝吸収線量とみなす
（README/[[lessons_learned]]の設計判断と同じ割り切り）。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

_MEV_TO_JOULE = 1.602176634e-13
_G_TO_KG = 1e-3


@dataclass
class VoxelGrid:
    origin_cm: np.ndarray          # (3,) グリッド原点（最小コーナー）
    shape: tuple                   # (nx, ny, nz)
    voxel_size_cm: float
    kerma_keV: np.ndarray = field(init=False)  # (nx,ny,nz) 積算カーマ [keV]

    def __post_init__(self):
        self.kerma_keV = np.zeros(self.shape, dtype=float)

    @classmethod
    def from_bbox(cls, bbox_min: np.ndarray, bbox_max: np.ndarray, resolution_cm: float) -> "VoxelGrid":
        extent = bbox_max - bbox_min
        shape = tuple(max(1, int(np.ceil(x / resolution_cm))) for x in extent)
        return cls(origin_cm=np.asarray(bbox_min, dtype=float), shape=shape, voxel_size_cm=resolution_cm)

    def voxel_index(self, points: np.ndarray):
        """点(N,3) -> ボクセル添字(N,3)とグリッド内かどうか(N,)。"""
        idx = np.floor((points - self.origin_cm) / self.voxel_size_cm).astype(int)
        shape_arr = np.array(self.shape)
        valid = np.all((idx >= 0) & (idx < shape_arr), axis=1)
        return idx, valid

    def voxel_centers(self) -> np.ndarray:
        """全ボクセル中心座標 (nx*ny*nz, 3)。"""
        nx, ny, nz = self.shape
        xs = self.origin_cm[0] + (np.arange(nx) + 0.5) * self.voxel_size_cm
        ys = self.origin_cm[1] + (np.arange(ny) + 0.5) * self.voxel_size_cm
        zs = self.origin_cm[2] + (np.arange(nz) + 0.5) * self.voxel_size_cm
        gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
        return np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)

    def voxel_volume_cm3(self) -> float:
        return self.voxel_size_cm ** 3

    def total_kerma_MeV(self) -> float:
        return float(self.kerma_keV.sum()) / 1000.0

    def dose_map_Gy(self, density_g_cm3: np.ndarray) -> np.ndarray:
        """材料密度マップ(shapeと同じ)からボクセルごとの吸収線量[Gy]を計算。"""
        mass_kg = density_g_cm3 * self.voxel_volume_cm3() * _G_TO_KG
        energy_J = self.kerma_keV * 1e-3 * _MEV_TO_JOULE
        return np.divide(energy_J, mass_kg, out=np.zeros_like(energy_J), where=mass_kg > 0)


def accumulate_track_length(grid: VoxelGrid, origin: np.ndarray, direction: np.ndarray,
                             length_cm: np.ndarray, weight_keV_per_cm: np.ndarray,
                             substep_cm: float | None = None, max_substeps: int = 40) -> None:
    """区間(origin, direction, length_cm)ごとに weight_keV_per_cm * dl をグリッドへ積算する。

    weight_keV_per_cm = E * (μen/ρ) * ρ は区間内で一定（材料・エネルギーとも不変の前提）。
    区間をサブステップに等分し、各サブステップの中点が属するボクセルへ
    weight_keV_per_cm * (length_cm/nsub) を加算する。
    """
    n = origin.shape[0]
    if n == 0:
        return
    if substep_cm is None:
        substep_cm = grid.voxel_size_cm / 2.0
    nsub = np.clip(np.ceil(length_cm / substep_cm).astype(int), 1, max_substeps)
    max_n = int(nsub.max())

    j = np.arange(max_n)
    frac = (j[None, :] + 0.5) / nsub[:, None]           # (n, max_n)
    valid = j[None, :] < nsub[:, None]                   # (n, max_n)

    points = (origin[:, None, :] + direction[:, None, :]
              * (length_cm[:, None] * frac)[:, :, None])  # (n, max_n, 3)
    sub_weight = weight_keV_per_cm * (length_cm / nsub)   # (n,)

    points_flat = points.reshape(-1, 3)
    weight_flat = np.broadcast_to(sub_weight[:, None], (n, max_n)).reshape(-1)
    valid_flat = valid.reshape(-1)

    idx, in_grid = grid.voxel_index(points_flat)
    keep = valid_flat & in_grid
    idx, weight_flat = idx[keep], weight_flat[keep]
    if len(idx) == 0:
        return
    flat_idx = np.ravel_multi_index((idx[:, 0], idx[:, 1], idx[:, 2]), grid.shape)
    np.add.at(grid.kerma_keV.reshape(-1), flat_idx, weight_flat)
