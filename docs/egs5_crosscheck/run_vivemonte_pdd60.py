"""viveMonte側の相互検証実行スクリプト（Phase 2b: 水ファントム深部線量分布PDD
＋2深さの側方プロファイル、60 keV単色）。

Phase 2a（BSF、docs/egs5_crosscheck/run_vivemonte_bsf60.py）と同一の手法
（TrajectoryRecorderの生の飛行区間から、スコア領域との解析的な厳密重なり長を
計算するtrack-length kerma推定量）を、単一の薄層ではなく複数の深さ・側方ビンに
拡張したもの。

幾何・物理条件（docs/plan_phits_crosscheck.md Phase 2の事前条件、Phase 2aと同一）:
- 60 keV単色光子、10×10cm²照射野、平行ビーム近似（半開角2.9°を無視）
- 水ファントム30×30×20cm、前面z=0
- 採点はファントム内部のみ（自由空気は比較対象にしない — Phase 2aの
  AIRDET密度トリック論争を踏まえ、PDDは絶対Gy/history同士の比較なので
  自由空気の参照値は不要）

ビン定義（EGS5側と厳密にビン境界を一致させる — 両コード共通の設計値）:
- PDD中心軸: |x|<=1cm, |y|<=1cm（断面2x2cm）、深さ0〜15cmを1cm刻み15ビン
- 側方プロファイル（表面近傍 z in [0,1]cm、10cm深 z in [9,10]cm）:
  |y|<=1cm、x方向を-8〜8cmの範囲で1cm刻み16ビン
  （表面近傍ビンはPDDの第1ビンと共通、10cm深ビンはPDDの第10ビンと共通）

出力は各ビンの吸収線量カーマ [Gy/history]（史ごとのΣx・Σx²から標準誤差を計算、
バッチ平均法ではない）。1 keV/g = 1.602176634e-13 Gy。
"""
from __future__ import annotations

import json

import numpy as np

from vivemonte.geometry import Geometry
from vivemonte.materials import density, mu_en_rho
from vivemonte.transport import transport_photons
from vivemonte.trajectory import TrajectoryRecorder

ENERGY_KEV = 60.0
FIELD_HALF_CM = 5.0          # 10x10 cm^2 照射野（平行ビーム近似）
N_TOTAL = 5_000_000
N_BATCH = 100_000
SEED = 1
KEV_PER_G_TO_GY = 1.602176634e-13

COL_HALF_CM = 1.0             # PDD中心軸コラムの半幅（|x|<=1, |y|<=1）
DEPTH_EDGES_CM = np.arange(0.0, 15.01, 1.0)   # 15ビン、1cm刻み
LATERAL_EDGES_CM = np.arange(-8.0, 8.01, 1.0)  # 16ビン、1cm刻み
LATERAL_DEPTHS = [("shallow", 0.0, 1.0), ("10cm", 9.0, 10.0)]


def _segment_box_overlap_cm(starts: np.ndarray, ends: np.ndarray,
                             lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """線分(starts[i]->ends[i])と軸並行直方体[lo,hi]の厳密な重なり長[cm]。

    docs/egs5_crosscheck/run_vivemonte_bsf60.py と同一実装（vive-auditor監査済み）。
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
    x = rng.uniform(-FIELD_HALF_CM, FIELD_HALF_CM, n)
    y = rng.uniform(-FIELD_HALF_CM, FIELD_HALF_CM, n)
    z = np.full(n, -1e-4)
    pos = np.column_stack([x, y, z])
    dirv = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    energy = np.full(n, ENERGY_KEV)
    return pos, dirv, energy


def _build_bins():
    """(名前, lo, hi) のリストを組み立てる。lo/hiはnp.array([x,y,z])。"""
    bins = []
    for i in range(len(DEPTH_EDGES_CM) - 1):
        z0, z1 = DEPTH_EDGES_CM[i], DEPTH_EDGES_CM[i + 1]
        bins.append((f"pdd_z{z0:.0f}-{z1:.0f}",
                     np.array([-COL_HALF_CM, -COL_HALF_CM, z0]),
                     np.array([COL_HALF_CM, COL_HALF_CM, z1])))
    for depth_name, z0, z1 in LATERAL_DEPTHS:
        for i in range(len(LATERAL_EDGES_CM) - 1):
            x0, x1 = LATERAL_EDGES_CM[i], LATERAL_EDGES_CM[i + 1]
            bins.append((f"lat_{depth_name}_x{x0:.0f}-{x1:.0f}",
                         np.array([x0, -COL_HALF_CM, z0]),
                         np.array([x1, COL_HALF_CM, z1])))
    return bins


def main() -> None:
    bins = _build_bins()
    print(f"エネルギー={ENERGY_KEV} keV, 照射野={2*FIELD_HALF_CM}x{2*FIELD_HALF_CM} cm^2 "
          f"(平行ビーム近似), n={N_TOTAL}, seed={SEED}, ビン数={len(bins)}")

    phantom_geom = Geometry([{
        "name": "phantom", "shape": "box", "material": "water",
        "center": [0.0, 0.0, 10.0], "size_cm": [30.0, 30.0, 20.0],
    }])

    rng = np.random.default_rng(SEED)
    per_history = {name: np.zeros(N_TOTAL) for name, _, _ in bins}
    done = 0
    while done < N_TOTAL:
        n = min(N_BATCH, N_TOTAL - done)
        pos, dirv, energy = _sample_parallel_beam(n, rng)
        recorder = TrajectoryRecorder()
        transport_photons(pos, dirv, energy, phantom_geom, rng, grid=None, recorder=recorder)

        starts = np.concatenate(recorder.starts)
        ends = np.concatenate(recorder.ends)
        seg_energy = np.concatenate(recorder.energies)
        photon_ids = np.concatenate(recorder.photon_ids)
        mu_en_linear = mu_en_rho("water", seg_energy) * density("water")

        for name, lo, hi in bins:
            overlap_cm = _segment_box_overlap_cm(starts, ends, lo, hi)
            hit = overlap_cm > 0
            if not np.any(hit):
                continue
            contrib_keV = overlap_cm[hit] * seg_energy[hit] * mu_en_linear[hit]
            batch_local = np.zeros(n)
            np.add.at(batch_local, photon_ids[hit], contrib_keV)
            per_history[name][done:done + n] = batch_local

        done += n
        print(f"  {done}/{N_TOTAL} histories 完了", flush=True)

    results = {}
    for name, lo, hi in bins:
        vol_cm3 = float(np.prod(hi - lo))
        mass_g = density("water") * vol_cm3
        arr_keV = per_history[name]
        mean_keV = arr_keV.mean()
        sem_keV = arr_keV.std(ddof=1) / np.sqrt(N_TOTAL)
        mean_gy = (mean_keV / mass_g) * KEV_PER_G_TO_GY
        sem_gy = (sem_keV / mass_g) * KEV_PER_G_TO_GY
        rel_err = sem_gy / mean_gy if mean_gy > 0 else float("nan")
        results[name] = {
            "lo_cm": lo.tolist(), "hi_cm": hi.tolist(), "mass_g": mass_g,
            "mean_Gy_per_history": mean_gy, "sem_Gy_per_history": sem_gy,
            "rel_err": rel_err,
        }
        print(f"[{name}] {mean_gy:.6e} +- {sem_gy:.6e} Gy/history "
              f"(相対{100*rel_err:.3f}%)")

    out_path = "docs/egs5_crosscheck/vivemonte_pdd60_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    print(f"\n結果を {out_path} に保存しました。")


if __name__ == "__main__":
    main()
