"""線量グリッドの後処理と「最大値が非物理的な位置に落ちていないか」の診断。

`chatcarlo run --dose-grid` の最大吸収線量・最大H*(10)は、背景（空気）ボクセルや
点線源の1/r²発散近傍に落ちることがある（docs/lessons_learned.md参照）。
CLIはここにある判定関数で該当ケースを検出して警告を出す。
"""
from __future__ import annotations

import numpy as np

from .geometry import Geometry
from .materials import density
from .tally import VoxelGrid


def dose_map_Gy(grid: VoxelGrid, geometry: Geometry) -> np.ndarray:
    """ボクセル中心の材料を判定し、その密度でカーマ→吸収線量[Gy]に換算する。

    グリッドはタリー専用であり材料を保持しないため、密度は出力時に
    ジオメトリーへ問い合わせて求める（ボクセル解像度が粗い場合、
    境界付近のボクセルは中心点1点で代表材料を決める近似になる）。
    """
    centers = grid.voxel_centers()
    mat = geometry.material_at(centers)
    density_flat = np.array([density(m) for m in mat])
    return grid.dose_map_Gy(density_flat.reshape(grid.shape))


def max_voxel_position_cm(grid: VoxelGrid, data: np.ndarray) -> np.ndarray:
    """dataの最大値を持つボクセルの中心座標[cm]。"""
    idx = np.unravel_index(int(np.argmax(data)), data.shape)
    return grid.origin_cm + (np.asarray(idx, dtype=float) + 0.5) * grid.voxel_size_cm


def background_medium_warning(material: str, background: str) -> str | None:
    """吸収線量の最大値ボクセルが背景（既定air）かどうかを判定する。

    吸収線量 = カーマ/密度 は媒質固有の量（同じカーマでも密度が違えば
    値が変わる）。空気は密度が非常に小さいため、材料境界のすぐ外側の
    空気ボクセルはカーマが同程度でも線量が大きく増幅されて見えることがある
    （[[lessons_learned]]参照）。この値は患者・検出器等、実体のある位置の
    被ばく評価には使えない。
    """
    if material != background:
        return None
    return ("最大値は空気中ボクセル（材料=air）です。吸収線量は媒質固有の量のため、"
            "この値は患者・検出器等の実体がある位置の被ばく評価には使えません。")


def near_source_air_warning(material: str, background: str, distance_from_source_cm: float,
                             nearest_object_distance_cm: float | None) -> str | None:
    """H*(10)最大値が点線源モデルの1/r²発散による非物理的な値かどうかを判定する。

    材料が背景（空気）かつ、シーン内のどの物体よりも線源に近い位置にある場合のみ
    警告する。「シーン内に実在する物体よりも線源に近い」は、その位置に人や検出器が
    存在し得ないことの明確な根拠になる（実際のX線管は housing/コリメータで
    覆われているが、ChatCarloの点線源モデルはそれを持たない）。
    """
    if material != background:
        return None
    if nearest_object_distance_cm is None or distance_from_source_cm >= nearest_object_distance_cm:
        return None
    return (f"最大値は線源から{distance_from_source_cm:.1f}cmの空気中ボクセルで、"
            f"シーン内のどの物体（最寄り{nearest_object_distance_cm:.1f}cm）よりも"
            "線源に近い位置です。点線源モデルの1/r²発散による非物理的な値であり、"
            "実在する位置の被ばく評価には使えません。評価したい位置（患者表面・"
            "操作者位置等）には直接細かいグリッドを敷いて計算してください。")
