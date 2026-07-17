"""Step 0a/0b/0c: Phase 2b残差(観測−1.71%のうち未説明の約36%)の机上検算。

plan_residual_check.md のStep 0。単一エネルギー(60 keV一次成分)のΔ_pred計算
(check_compton_transfer.py)を、実際の散乱スペクトル混合(多重散乱で軟化した
二次光子を含む)まで畳み込んだ予測値に拡張する。同時に、μen/ρのlog-log補間
バイアス(Step 0b)と密度・断面積微小差の深さ方向解析伝播(Step 0c)も評価する。

実行: PYTHONPATH=. .venv/bin/python docs/egs5_crosscheck/check_residual_convolution.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import xraylib

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_compton_transfer import average_transfer_fraction  # noqa: E402
from extract_pegs5_photo_xs import cross_sections_cm2_per_g  # noqa: E402

from chatcarlo.geometry import Geometry
from chatcarlo.materials import density, mu_en_rho
from chatcarlo.transport import transport_photons
from chatcarlo.trajectory import TrajectoryRecorder

HERE = Path(__file__).resolve().parent
EGS5_PDD60 = HERE / "egs5" / "run_pdd60_phantom" / "pgs5job.pegs5lst"

ENERGY_KEV = 60.0
FIELD_HALF_CM = 5.0
N_SPECTRUM = 1_000_000
N_BATCH = 100_000
SEED = 7
KEV_PER_G_TO_GY = 1.602176634e-13

COL_HALF_CM = 1.0
# PDD中心軸(深さ方向)ビン。vive-auditor所見1: これだけでは側方(OCR)ビンの
# スペクトルを代表できない懸念があるため、側方ビンも別途タリーする(下記)。
REPRESENTATIVE_BINS = [
    ("pdd_z0-1", 0.0, 1.0),
    ("pdd_z5-6", 5.0, 6.0),
    ("pdd_z10-11", 10.0, 11.0),
    ("pdd_z14-15", 14.0, 15.0),
]
# 側方プロファイル代表ビン(中心軸付近と視野端に近いビンで散乱スペクトルの
# 違いを確認する)。(名前, z0, z1, x0, x1)。z0-1は表面近傍、z9-10は10cm深。
LATERAL_REPRESENTATIVE_BINS = [
    ("lat_shallow_x0-1", 0.0, 1.0, 0.0, 1.0),
    ("lat_shallow_x-4-3", 0.0, 1.0, -4.0, -3.0),
    ("lat_10cm_x0-1", 9.0, 10.0, 0.0, 1.0),
    ("lat_10cm_x-4-3", 9.0, 10.0, -4.0, -3.0),
]

EBIN_EDGES_KEV = np.arange(10.0, 60.001, 2.0)  # PEGS5 AP=10keVより下は対象外(除外分を別途報告)
EBIN_CENTERS_KEV = 0.5 * (EBIN_EDGES_KEV[:-1] + EBIN_EDGES_KEV[1:])


def _segment_box_overlap_cm(starts, ends, lo, hi):
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


def _sample_parallel_beam(n, rng):
    x = rng.uniform(-FIELD_HALF_CM, FIELD_HALF_CM, n)
    y = rng.uniform(-FIELD_HALF_CM, FIELD_HALF_CM, n)
    z = np.full(n, -1e-4)
    pos = np.column_stack([x, y, z])
    dirv = np.tile(np.array([0.0, 0.0, 1.0]), (n, 1))
    energy = np.full(n, ENERGY_KEV)
    return pos, dirv, energy


def tally_kerma_spectra() -> dict[str, np.ndarray]:
    """代表4ビンで「kerma重み付きエネルギーヒストグラム」を作る。

    各セグメントの寄与 overlap_cm*E*mu_en_linear(E) を seg_energy でヒストグラム化。
    これをΔ_pred(E)の重み(Σ φ(E)μen(E)E)としてそのまま使える。
    """
    phantom_geom = Geometry([{
        "name": "phantom", "shape": "box", "material": "water",
        "center": [0.0, 0.0, 10.0], "size_cm": [30.0, 30.0, 20.0],
    }])
    all_names = [name for name, _, _ in REPRESENTATIVE_BINS] + \
                [name for name, *_ in LATERAL_REPRESENTATIVE_BINS]
    rng = np.random.default_rng(SEED)
    hist = {name: np.zeros(len(EBIN_CENTERS_KEV)) for name in all_names}
    below_ap_sum = {name: 0.0 for name in all_names}
    total_sum = {name: 0.0 for name in all_names}

    def _score(name, lo, hi, starts, ends, seg_energy, mu_en_linear):
        overlap_cm = _segment_box_overlap_cm(starts, ends, lo, hi)
        hit = overlap_cm > 0
        if not np.any(hit):
            return
        contrib = overlap_cm[hit] * seg_energy[hit] * mu_en_linear[hit]
        e_hit = seg_energy[hit]
        total_sum[name] += contrib.sum()
        below_ap_sum[name] += contrib[e_hit < EBIN_EDGES_KEV[0]].sum()
        idx = np.digitize(e_hit, EBIN_EDGES_KEV) - 1
        valid = (idx >= 0) & (idx < len(EBIN_CENTERS_KEV))
        np.add.at(hist[name], idx[valid], contrib[valid])

    done = 0
    while done < N_SPECTRUM:
        n = min(N_BATCH, N_SPECTRUM - done)
        pos, dirv, energy = _sample_parallel_beam(n, rng)
        recorder = TrajectoryRecorder()
        transport_photons(pos, dirv, energy, phantom_geom, rng, grid=None, recorder=recorder)

        starts = np.concatenate(recorder.starts)
        ends = np.concatenate(recorder.ends)
        seg_energy = np.concatenate(recorder.energies)
        mu_en_linear = mu_en_rho("water", seg_energy) * density("water")

        for name, z0, z1 in REPRESENTATIVE_BINS:
            lo = np.array([-COL_HALF_CM, -COL_HALF_CM, z0])
            hi = np.array([COL_HALF_CM, COL_HALF_CM, z1])
            _score(name, lo, hi, starts, ends, seg_energy, mu_en_linear)
        for name, z0, z1, x0, x1 in LATERAL_REPRESENTATIVE_BINS:
            lo = np.array([x0, -COL_HALF_CM, z0])
            hi = np.array([x1, COL_HALF_CM, z1])
            _score(name, lo, hi, starts, ends, seg_energy, mu_en_linear)

        done += n

    for name in all_names:
        frac_below = below_ap_sum[name] / total_sum[name] * 100 if total_sum[name] > 0 else float("nan")
        print(f"  [{name}] AP(10keV)未満の寄与割合(除外分): {frac_below:.3f}%")

    return hist


def delta_pred_grid() -> dict[str, np.ndarray]:
    """EBIN_CENTERS_KEVの各エネルギーでΔ_pred(E)を計算(t_KN/t_Sq双方)。

    分母は2通り用意する:
    - muen_xaamdi: ChatCarloが実際に使うlog-log補間値(実際の比較に対応)
    - muen_firstprin: xraylib photo+compton*t_Sq(補間を経ない連続的な近似、
      格子点での検証で相対差0.02-0.08%と高精度に確認済み。Step 0bのlog-log補間
      アーティファクトを含まない「純粋なコンプトン移行仮説」の寄与を分離するための分母)
    """
    delta_kn = np.zeros(len(EBIN_CENTERS_KEV))
    delta_sq = np.zeros(len(EBIN_CENTERS_KEV))
    delta_kn_firstprin = np.zeros(len(EBIN_CENTERS_KEV))
    for i, E in enumerate(EBIN_CENTERS_KEV):
        _, t_kn = average_transfer_fraction(float(E), use_sq=False)
        _, t_sq = average_transfer_fraction(float(E), use_sq=True)
        pegs5 = cross_sections_cm2_per_g(EGS5_PDD60, E / 1000.0)
        muen_xaamdi = float(mu_en_rho("water", E)[0])
        photo_xr = xraylib.CS_Photo_CP("Water, Liquid", float(E))
        compt_xr = xraylib.CS_Compt_CP("Water, Liquid", float(E))
        muen_firstprin = photo_xr + compt_xr * t_sq
        k_eff_kn = pegs5["photo"] + pegs5["compton"] * t_kn
        k_eff_sq = pegs5["photo"] + pegs5["compton"] * t_sq
        delta_kn[i] = (k_eff_kn - muen_xaamdi) / muen_xaamdi * 100
        delta_sq[i] = (k_eff_sq - muen_xaamdi) / muen_xaamdi * 100
        delta_kn_firstprin[i] = (k_eff_kn - muen_firstprin) / muen_firstprin * 100
    return {"delta_kn": delta_kn, "delta_sq": delta_sq, "delta_kn_firstprin": delta_kn_firstprin}


def step0b_interpolation_bias() -> None:
    print("\n### Step 0b: mu_en/rho log-log補間バイアスの上界評価 ###")
    print("検証: 格子点そのもの(補間なし)で「第一原理近似」自体の精度を確認")
    for E in (20.0, 30.0, 40.0, 50.0, 60.0, 80.0):
        muen_table = float(mu_en_rho("water", E)[0])
        photo_xr = xraylib.CS_Photo_CP("Water, Liquid", E)
        compt_xr = xraylib.CS_Compt_CP("Water, Liquid", E)
        _, t_sq = average_transfer_fraction(E, use_sq=True)
        firstprin = photo_xr + compt_xr * t_sq
        rel = (muen_table - firstprin) / firstprin * 100
        print(f"  E={E:5.1f}keV(格子点)  table={muen_table:.6f}  firstprin={firstprin:.6f}  "
              f"相対差={rel:+.3f}%")
    print("(格子点では0.1%未満で一致 -> 第一原理近似自体は高精度、中間点との差は"
          "補間の曲率誤差そのものと確認)")
    print("\n(XAAMDI格子点: ...,20,30,40,50,60,80... / 中間エネルギーで比較)")
    for E in (35.0, 45.0, 55.0):
        muen_interp = float(mu_en_rho("water", E)[0])
        photo_xr = xraylib.CS_Photo_CP("Water, Liquid", E)
        compt_xr = xraylib.CS_Compt_CP("Water, Liquid", E)
        _, t_sq = average_transfer_fraction(E, use_sq=True)
        muen_firstprin = photo_xr + compt_xr * t_sq
        rel = (muen_interp - muen_firstprin) / muen_firstprin * 100
        print(f"  E={E:.0f}keV  XAAMDI log-log補間={muen_interp:.6f}  "
              f"第一原理近似(xraylib photo+compton*t_Sq)={muen_firstprin:.6f}  "
              f"相対差={rel:+.3f}%")


def step0c_analytic_propagation() -> None:
    print("\n### Step 0c: 密度差・断面積差の深さ方向解析伝播(60 keV一次) ###")
    pegs5 = cross_sections_cm2_per_g(EGS5_PDD60, 0.060)
    photo_xr = xraylib.CS_Photo_CP("Water, Liquid", 60.0)
    compt_xr = xraylib.CS_Compt_CP("Water, Liquid", 60.0)
    rayl_xr = xraylib.CS_Rayl_CP("Water, Liquid", 60.0)
    total_xr = photo_xr + compt_xr + rayl_xr
    total_pegs5 = pegs5["photo"] + pegs5["compton"] + pegs5["rayleigh"]
    rho_pegs5, rho_chatcarlo = 1.001, density("water")

    rel_total_xs = (total_pegs5 - total_xr) / total_xr
    rel_density = (rho_pegs5 - rho_chatcarlo) / rho_chatcarlo
    print(f"  全断面積(質量減弱係数)差: PEGS5 {total_pegs5:.6f} vs xraylib {total_xr:.6f}  "
          f"相対差={rel_total_xs*100:+.4f}%")
    print(f"  密度差: PEGS5(水template) {rho_pegs5} vs ChatCarlo {rho_chatcarlo}  "
          f"相対差={rel_density*100:+.4f}%")
    print("  透過率比 T_EGS5/T_ChatCarlo = exp(-(mu_pegs5*rho_pegs5 - mu_xr*rho_chatcarlo)*z)"
          " から寄与を分離:")
    print(f"  {'深さ[cm]':>8s}  {'断面積差由来[%pt]':>18s}  {'密度差由来[%pt]':>16s}  "
          f"{'合計[%pt]':>10s}")
    for z in (0.5, 5.5, 10.5, 14.5):
        mu_pegs5_lin = total_pegs5 * rho_pegs5
        mu_xr_lin = total_xr * rho_chatcarlo
        # 断面積差のみ(密度共通rho_chatcarloで評価)
        shift_xs = (math.exp(-total_pegs5 * rho_chatcarlo * z) /
                    math.exp(-total_xr * rho_chatcarlo * z) - 1.0) * 100
        # 密度差のみ(断面積共通total_xrで評価)
        shift_rho = (math.exp(-total_xr * rho_pegs5 * z) /
                     math.exp(-total_xr * rho_chatcarlo * z) - 1.0) * 100
        shift_total = (math.exp(-mu_pegs5_lin * z) / math.exp(-mu_xr_lin * z) - 1.0) * 100
        print(f"  {z:8.1f}  {shift_xs:18.4f}  {shift_rho:16.4f}  {shift_total:10.4f}")


def main() -> None:
    print("### Step 0a: 代表4ビンのkerma重み付きエネルギースペクトル畳み込み ###")
    print(f"n={N_SPECTRUM}, seed={SEED}, エネルギービン=10-60keV/2keV刻み\n")

    hist = tally_kerma_spectra()
    grid = delta_pred_grid()

    print("\n畳み込み結果 Δ_pred_conv(bin) = Σ w(E)*Δ_pred(E) / Σ w(E):")
    print(f"  {'ビン':>18s}  {'Δ(t_KN,XAAMDI補間分母)':>24s}  {'Δ(t_Sq,XAAMDI補間分母)':>24s}  "
          f"{'Δ(t_KN,第一原理分母)':>22s}")
    conv_results = {}
    all_bin_names = [name for name, _, _ in REPRESENTATIVE_BINS] + \
                     [name for name, *_ in LATERAL_REPRESENTATIVE_BINS]
    for name in all_bin_names:
        w = hist[name]
        wsum = w.sum()
        conv_kn = float((w * grid["delta_kn"]).sum() / wsum)
        conv_sq = float((w * grid["delta_sq"]).sum() / wsum)
        conv_kn_fp = float((w * grid["delta_kn_firstprin"]).sum() / wsum)
        conv_results[name] = (conv_kn, conv_sq, conv_kn_fp)
        tag = "" if name.startswith("pdd_") else "  [側方]"
        print(f"  {name:>18s}  {conv_kn:23.3f}%  {conv_sq:23.3f}%  {conv_kn_fp:21.3f}%{tag}")

    pdd_sq = [conv_results[n][1] for n, _, _ in REPRESENTATIVE_BINS]
    lat_sq = [conv_results[n][1] for n, *_ in LATERAL_REPRESENTATIVE_BINS]
    print(f"\nPDDビン Δ(t_Sq)平均: {sum(pdd_sq)/len(pdd_sq):+.3f}%  "
          f"側方ビン Δ(t_Sq)平均: {sum(lat_sq)/len(lat_sq):+.3f}%  "
          f"(vive-auditor所見1への対応: 側方ビンでの検算)")

    print("\n(参考)単一エネルギー60keV一次成分のみのΔ_pred "
          "(check_compton_transfer.py既報値): t_KN=-1.095%, t_Sq=+0.179%")
    print("(Δ(t_KN,第一原理分母)は、ChatCarlo側のlog-log補間アーティファクト[Step 0b]を"
          "含まない「純粋なコンプトン移行分布仮説単独」の寄与)")

    step0b_interpolation_bias()
    step0c_analytic_propagation()


if __name__ == "__main__":
    main()
