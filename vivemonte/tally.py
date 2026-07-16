"""ボクセル吸収線量・周辺線量当量タリー — 輸送ジオメトリーとは独立なグリッド。

輸送は解析面トラッキング（vivemonte/geometry.py, vivemonte/transport.py）
で行い、スコアリングだけをこのモジュールの一様グリッドに刻む。

track-length estimator（trackごとの飛程積分）で2種類の量を積算する:
  カーマ:    K += E * (μen/ρ) * ρ * dl                      [keV]
  H*(10):   H += (h*(10)/Φ)(E) * dl、後でボクセル体積で正規化   [pSv]
どちらも区間内はエネルギー・材料とも一定なので積分自体に離散化誤差はない。
「その区間がどのボクセルに何cm分入っているか」の空間分配は層化乱数点による
モンテカルロ分配で、任意のサブステップ長で不偏（期待値＝厳密な重なり長。
サブステップ長は分散にのみ影響する。accumulate_track_lengthのdocstring参照）。
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
    kerma_keV: np.ndarray = field(init=False)      # (nx,ny,nz) 積算カーマ [keV]
    h10_track_pSv_cm3: np.ndarray = field(init=False)  # (nx,ny,nz) H*(10)飛程積分 [pSv・cm³]（体積正規化前）

    def __post_init__(self):
        self.kerma_keV = np.zeros(self.shape, dtype=float)
        self.h10_track_pSv_cm3 = np.zeros(self.shape, dtype=float)

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

    def h10_map_pSv(self) -> np.ndarray:
        """ボクセルごとの周辺線量当量H*(10) [pSv]（飛程積分をボクセル体積で正規化）。"""
        return self.h10_track_pSv_cm3 / self.voxel_volume_cm3()


def accumulate_track_length(target: np.ndarray, grid: VoxelGrid, origin: np.ndarray,
                             direction: np.ndarray, length_cm: np.ndarray, weight_per_cm: np.ndarray,
                             rng: np.random.Generator,
                             substep_cm: float | None = None, max_substeps: int = 40) -> None:
    """区間(origin, direction, length_cm)ごとに weight_per_cm * dl を target グリッドへ積算する。

    target は grid.shape と同じ形の任意の量（カーマ・H*(10)飛程積分など）で、
    weight_per_cm は区間内で一定（材料・エネルギーとも不変の前提）とする。
    区間をサブステップに等分し、各サブステップ内の一様乱数点（層化サンプリング）が
    属するボクセルへ weight_per_cm * (length_cm/nsub) を加算する。

    乱数点を使う理由: 以前はサブステップの中点（決定的）を使っていたが、
    区間の始点がボクセル境界ちょうどに揃う条件（例: parallel照射野で全光子が
    ファントム前面から出発）では量子化誤差の位相が全区間で同期し、表面ボクセルで
    約-3%の系統的過小評価になる（独立監査で発見）。層化乱数点なら任意のボクセルに
    対して期待値が厳密な幾何学的重なり長に一致する（不偏推定量）。
    rng は輸送用と別のストリームを渡すこと（タリーが輸送の乱数列を消費して
    物理結果を変えないため — transport_photons が spawn で自動生成する）。
    """
    n = origin.shape[0]
    if n == 0:
        return
    if substep_cm is None:
        substep_cm = grid.voxel_size_cm / 2.0
    nsub = np.clip(np.ceil(length_cm / substep_cm).astype(int), 1, max_substeps)
    max_n = int(nsub.max())

    j = np.arange(max_n)
    frac = (j[None, :] + rng.random((n, max_n))) / nsub[:, None]  # (n, max_n) 層化乱数点
    valid = j[None, :] < nsub[:, None]                   # (n, max_n)

    points = (origin[:, None, :] + direction[:, None, :]
              * (length_cm[:, None] * frac)[:, :, None])  # (n, max_n, 3)
    sub_weight = weight_per_cm * (length_cm / nsub)       # (n,)

    points_flat = points.reshape(-1, 3)
    weight_flat = np.broadcast_to(sub_weight[:, None], (n, max_n)).reshape(-1)
    valid_flat = valid.reshape(-1)

    idx, in_grid = grid.voxel_index(points_flat)
    keep = valid_flat & in_grid
    idx, weight_flat = idx[keep], weight_flat[keep]
    if len(idx) == 0:
        return
    flat_idx = np.ravel_multi_index((idx[:, 0], idx[:, 1], idx[:, 2]), grid.shape)
    np.add.at(target.reshape(-1), flat_idx, weight_flat)
