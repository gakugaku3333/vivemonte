"""viveMonte側の相互検証実行スクリプト（Phase 2a: 後方散乱係数BSF、60 keV単色）。

EGS5側（docs/egs5_crosscheck/bsf60_free/, bsf60_phantom/、詳細は
docs/egs5_crosscheck/bsf60_NOTES.md）と同一幾何・物理条件でBSFを計算する:

- 60 keV単色光子、10×10cm²照射野（平行ビーム近似、EGS5側と同じ簡略化）
- 自由空気カーマ（分母）: 空気のみ、深さ0〜0.2cm層でスコア
- 水ファントム表面層カーマ（分子、後方散乱込み）: 水30×30×20cm、
  前面から深さ0〜0.2cm層でスコア
- BSF = 表面層カーマ ÷ 自由空気カーマ

実行:
    PYTHONPATH=. .venv/bin/python docs/egs5_crosscheck/run_vivemonte_bsf60.py

## なぜ既定の --dose-grid（VoxelGrid）を使わないか

VoxelGrid.accumulate_track_length は飛行区間をサブステップ分割して積算するが、
分割数は `substep_cm` あたり最大40分割に打ち切られる（vivemonte/tally.py）。
今回のスコア層は厚さ0.2cmしかないのに対し、水60 keVの平均自由行程は約5cmで、
1回目の飛行区間の多くはその何倍も長い。長い区間を高々40分割すると、0.2cm厚の
層に落ちる分割点がわずか1〜2個になり、真の重なり長（連続値）が粗い整数個の
分割単位に丸められる（CLAUDE.mdの「解像度を細かくするほど値が増大する」という
既知の境界効果と同根の問題）。

そこで本スクリプトでは TrajectoryRecorder が記録する生の飛行区間
（開始点・終了点・区間内一定のエネルギー）から、スコア領域（軸並行の直方体）
との厳密な解析的重なり長を計算する（区間ごとに1回のtの範囲計算、離散化なし）。
これは薄層カーマのtrack-length推定量として離散化誤差を持たない。副産物として
photon_idごとに正確な史あたりカーマが得られるため、EGS5と同じ史ごとの
Σx・Σx²モーメント統計で標準誤差を計算できる（バッチ平均法のような近似は不要）。

自由空気側の密度は実際の空気密度をそのまま使う（EGS5側が薄層での相互作用希少性
のためAIRDET密度×1000の分散低減トリックを要したのとは対照的——track-length
推定量は「相互作用が起きたか」ではなく解析的な区間重なり長で決まるため、密度を
上げて相互作用頻度を稼ぐ必要がない）。
"""
from __future__ import annotations

import numpy as np

from vivemonte.geometry import Geometry
from vivemonte.materials import density, mu_en_rho
from vivemonte.transport import transport_photons
from vivemonte.trajectory import TrajectoryRecorder

ENERGY_KEV = 60.0
FIELD_HALF_CM = 5.0          # 10x10 cm^2 照射野（平行ビーム近似）
LAYER_CM = 0.2               # スコア層厚（EGS5側 zbound/zlayer と同一）
N_TOTAL = 500_000            # EGS5側と同数（直接比較のため）
N_BATCH = 25_000
SEED = 1

# スコア領域（両ランで共通、軸並行の直方体）: |x|<=5, |y|<=5, 0<=z<=0.2
_LO = np.array([-FIELD_HALF_CM, -FIELD_HALF_CM, 0.0])
_HI = np.array([FIELD_HALF_CM, FIELD_HALF_CM, LAYER_CM])


def _segment_box_overlap_cm(starts: np.ndarray, ends: np.ndarray,
                             lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """線分(starts[i]->ends[i])と軸並行直方体[lo,hi]の厳密な重なり長[cm]。

    slab法（geometry.py の _intersect_box と同じ原理）をt∈[0,1]（線分自身の
    範囲）に制限して適用する。サブステップ分割を一切使わない解析解。
    """
    d = ends - starts
    length = np.linalg.norm(d, axis=1)
    t_enter = np.zeros(len(starts))
    t_exit = np.ones(len(starts))
    for k in range(3):
        dk, sk = d[:, k], starts[:, k]
        parallel = np.abs(dk) < 1e-12
        dk_safe = np.where(parallel, 1.0, dk)
        ta = (lo[k] - sk) / dk_safe
        tb = (hi[k] - sk) / dk_safe
        tmin_k = np.where(parallel, -np.inf, np.minimum(ta, tb))
        tmax_k = np.where(parallel, np.inf, np.maximum(ta, tb))
        outside_slab = parallel & ((sk < lo[k]) | (sk > hi[k]))
        tmin_k = np.where(outside_slab, 2.0, tmin_k)
        tmax_k = np.where(outside_slab, -1.0, tmax_k)
        t_enter = np.maximum(t_enter, tmin_k)
        t_exit = np.minimum(t_exit, tmax_k)
    overlap_t = np.clip(t_exit - t_enter, 0.0, None)
    return overlap_t * length


def _sample_parallel_beam(n: int, rng: np.random.Generator):
    """10x10cm^2照射野、平行ビーム（非発散近似）、60 keV単色。"""
    x = rng.uniform(-FIELD_HALF_CM, FIELD_HALF_CM, n)
    y = rng.uniform(-FIELD_HALF_CM, FIELD_HALF_CM, n)
    z = np.full(n, -1e-4)  # 前面のごくわずか手前から入射（境界上ちょうどを避ける）
    pos = np.column_stack([x, y, z])
    dirv = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    energy = np.full(n, ENERGY_KEV)
    return pos, dirv, energy


def run_side(name: str, geometry: Geometry, material: str, seed: int):
    """自由空気 or ファントムのどちらかのランを実行し、史ごとのカーマ配列を返す。"""
    rng = np.random.default_rng(seed)
    per_history = np.zeros(N_TOTAL)
    done = 0
    while done < N_TOTAL:
        n = min(N_BATCH, N_TOTAL - done)
        pos, dirv, energy = _sample_parallel_beam(n, rng)
        recorder = TrajectoryRecorder()
        transport_photons(pos, dirv, energy, geometry, rng, grid=None, recorder=recorder)

        starts = np.concatenate(recorder.starts)
        ends = np.concatenate(recorder.ends)
        seg_energy = np.concatenate(recorder.energies)
        photon_ids = np.concatenate(recorder.photon_ids)

        overlap_cm = _segment_box_overlap_cm(starts, ends, _LO, _HI)
        hit = overlap_cm > 0
        if np.any(hit):
            mu_en_linear = mu_en_rho(material, seg_energy[hit]) * density(material)
            contrib_keV = overlap_cm[hit] * seg_energy[hit] * mu_en_linear
            batch_local = np.zeros(n)
            np.add.at(batch_local, photon_ids[hit], contrib_keV)
            per_history[done:done + n] = batch_local

        done += n

    mean_keV = per_history.mean()
    sem_keV = per_history.std(ddof=1) / np.sqrt(N_TOTAL)
    print(f"[{name}] mean={mean_keV:.6e} keV/history  "
          f"SEM={sem_keV:.6e} keV/history  "
          f"relSEM={100 * sem_keV / mean_keV:.4f}%")
    return per_history, mean_keV, sem_keV


def main() -> None:
    print(f"エネルギー={ENERGY_KEV} keV, 照射野={2*FIELD_HALF_CM}x{2*FIELD_HALF_CM} cm^2 "
          f"(平行ビーム近似), スコア層厚={LAYER_CM} cm, n={N_TOTAL}/ラン, seed={SEED}")
    print(f"density(water)={density('water'):.4f} g/cm^3  "
          f"density(air)={density('air'):.6f} g/cm^3")
    print(f"mu_en/rho(water,60keV)={float(mu_en_rho('water', [ENERGY_KEV])[0]):.5f} cm^2/g  "
          f"mu_en/rho(air,60keV)={float(mu_en_rho('air', [ENERGY_KEV])[0]):.5f} cm^2/g")
    print()

    # 自由空気ジオメトリー: 実際の空気密度のまま（分散低減トリック不要、上記docstring参照）
    air_geom = Geometry([{
        "name": "air_slab", "shape": "box", "material": "air",
        "center": [0.0, 0.0, 2.5], "size_cm": [20.0, 20.0, 5.0],
    }])
    # ファントムジオメトリー: 水30x30x20cm、前面z=0（EGS5側と同一寸法）
    phantom_geom = Geometry([{
        "name": "phantom", "shape": "box", "material": "water",
        "center": [0.0, 0.0, 10.0], "size_cm": [30.0, 30.0, 20.0],
    }])

    free_hist, free_mean, free_sem = run_side("自由空気(分母)", air_geom, "air", SEED)
    phantom_hist, phantom_mean, phantom_sem = run_side("ファントム表面層(分子)", phantom_geom, "water", SEED + 1)

    # 質量あたりカーマへ換算（スコア層の質量=density*volume、EGS5側NOTES.mdと同じ手順）。
    # 生のエネルギー和のままでは空気層と水層で質量が約830倍違うため比較にならない。
    volume_cm3 = (2 * FIELD_HALF_CM) * (2 * FIELD_HALF_CM) * LAYER_CM
    mass_air_g = density("air") * volume_cm3
    mass_water_g = density("water") * volume_cm3
    kerma_air = free_mean / mass_air_g
    kerma_water = phantom_mean / mass_water_g
    print()
    print(f"スコア層体積={volume_cm3:.1f} cm^3  質量(air)={mass_air_g:.4f} g  質量(water)={mass_water_g:.4f} g")
    print(f"質量あたりカーマ: 自由空気={kerma_air:.6e} keV/g/history  "
          f"水表面層={kerma_water:.6e} keV/g/history")

    bsf = kerma_water / kerma_air
    rel_free = free_sem / free_mean
    rel_phantom = phantom_sem / phantom_mean
    bsf_rel_err = np.sqrt(rel_free ** 2 + rel_phantom ** 2)
    bsf_err = bsf * bsf_rel_err

    print()
    print(f"BSF (viveMonte) = {bsf:.4f} ± {bsf_err:.4f}  (1sigma, 相対{100*bsf_rel_err:.3f}%)")
    print()
    print("EGS5側 (docs/egs5_crosscheck/bsf60_NOTES.md, IBOUND=1, n=500,000/run):")
    print("  BSF (EGS5) = 1.4529 ± 0.0243  (1sigma, 相対1.670%)")
    egs5_bsf, egs5_err = 1.4529, 0.0243
    diff = bsf - egs5_bsf
    combined_err = np.sqrt(bsf_err ** 2 + egs5_err ** 2)
    n_sigma = abs(diff) / combined_err if combined_err > 0 else float("nan")
    print(f"  差 = {diff:+.4f}  (相対{100*diff/egs5_bsf:+.2f}%,  "
          f"合成誤差{combined_err:.4f}に対し{n_sigma:.2f}sigma)")


if __name__ == "__main__":
    main()
