"""ChatCarlo側の相互検証実行スクリプト（Phase 2a Step 2: BSF_w、60 keV単色）。

docs/egs5_crosscheck/plan_bsf60_resolution.md のStep 2に対応する。従来定義の
BSF（自由空気カーマを分母とする）は、EGS5側の密度トリック（AIRDET密度スケーリング）
に約3σの未解明な系統性が見つかった（BSF60_RESULTS.md「Step 1」参照）ため、
密度トリック自体を不要にする分母の再定義を使う:

    BSF_w = （水30x30x20cmファントム前面0〜0.2cm層のカーマ、後方散乱込み）
            ÷ （同一位置に置いた水30x30x0.2cm薄層“のみ”のカーマ、後方バルクなし）

分子・分母とも実密度の水（分散低減トリック不要）。差分は純粋に「後方バルクからの
後方散乱」であり、BSFが測りたい物理そのものが残る。BSF_wは文献の従来定義BSFとは
値が異なる（水/空気のμen/ρ比の分だけずれる）ので、文献値と直接比較しないこと
（文献照合はStep 3で別途、従来定義BSFで実施する）。

実行:
    PYTHONPATH=. .venv/bin/python docs/egs5_crosscheck/run_chatcarlo_bsf60_thinslab.py
"""
from __future__ import annotations

import numpy as np

from chatcarlo.geometry import Geometry
from chatcarlo.materials import density, mu_en_rho
from chatcarlo.transport import transport_photons
from chatcarlo.trajectory import TrajectoryRecorder

ENERGY_KEV = 60.0
FIELD_HALF_CM = 5.0          # 10x10 cm^2 照射野（平行ビーム近似）
LAYER_CM = 0.2               # スコア層厚（EGS5側と同一）
N_TOTAL = 500_000
N_BATCH = 25_000
SEED = 1

_LO = np.array([-FIELD_HALF_CM, -FIELD_HALF_CM, 0.0])
_HI = np.array([FIELD_HALF_CM, FIELD_HALF_CM, LAYER_CM])


def _segment_box_overlap_cm(starts: np.ndarray, ends: np.ndarray,
                             lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """線分(starts[i]->ends[i])と軸並行直方体[lo,hi]の厳密な重なり長[cm]。"""
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
    z = np.full(n, -1e-4)
    pos = np.column_stack([x, y, z])
    dirv = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    energy = np.full(n, ENERGY_KEV)
    return pos, dirv, energy


def run_side(name: str, geometry: Geometry, material: str, seed: int):
    """薄層のみ or フルファントムのどちらかのランを実行し、史ごとのカーマ配列を返す。"""
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
    print(f"density(water)={density('water'):.4f} g/cm^3")
    print()

    # 薄層のみ（後方バルクなし）: 水30x30x0.2cm、前面z=0（分母）
    thinslab_geom = Geometry([{
        "name": "thinslab", "shape": "box", "material": "water",
        "center": [0.0, 0.0, LAYER_CM / 2], "size_cm": [30.0, 30.0, LAYER_CM],
    }])
    # フルファントム: 水30x30x20cm、前面z=0（EGS5側と同一寸法、分子）
    phantom_geom = Geometry([{
        "name": "phantom", "shape": "box", "material": "water",
        "center": [0.0, 0.0, 10.0], "size_cm": [30.0, 30.0, 20.0],
    }])

    thin_hist, thin_mean, thin_sem = run_side("水薄層のみ(分母)", thinslab_geom, "water", SEED)
    phantom_hist, phantom_mean, phantom_sem = run_side("ファントム表面層(分子)", phantom_geom, "water", SEED + 1)

    volume_cm3 = (2 * FIELD_HALF_CM) * (2 * FIELD_HALF_CM) * LAYER_CM
    mass_water_g = density("water") * volume_cm3
    kerma_thin = thin_mean / mass_water_g
    kerma_phantom = phantom_mean / mass_water_g
    print()
    print(f"スコア層体積={volume_cm3:.1f} cm^3  質量(water)={mass_water_g:.4f} g")
    print(f"質量あたりカーマ: 水薄層のみ={kerma_thin:.6e} keV/g/history  "
          f"ファントム表面層={kerma_phantom:.6e} keV/g/history")

    bsf_w = kerma_phantom / kerma_thin
    rel_thin = thin_sem / thin_mean
    rel_phantom = phantom_sem / phantom_mean
    bsf_rel_err = np.sqrt(rel_thin ** 2 + rel_phantom ** 2)
    bsf_err = bsf_w * bsf_rel_err

    print()
    print(f"BSF_w (ChatCarlo) = {bsf_w:.4f} ± {bsf_err:.4f}  (1sigma, 相対{100*bsf_rel_err:.3f}%)")


if __name__ == "__main__":
    main()
