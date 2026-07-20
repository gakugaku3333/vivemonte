"""光子輸送カーネル — 解析面トラッキング（Woodcock不要）。

各光子は「次に消費すべき光学的厚み τ = -ln(ξ)」を持ち、材料境界を
跨ぐたびに τ から μ·Δs を差し引いていく。τ が現在の区間内で尽きたら
そこが実際の相互作用点になる。均質な区間ごとに μ が一定なので
Woodcock delta-trackingの仮想衝突は不要（空気の広い空間で無駄がない）。

周辺の責務は分離してある:
- スペクトル生成（SpekPy/Kramers・ヒール軸外スペクトル）: chatcarlo/spectrum.py
- 線源サンプリング・mAs光子数校正: chatcarlo/source.py
- 相互作用の角度・エネルギー抽選（束縛コンプトン/レイリー）: chatcarlo/physics.py
- 軌跡記録（trace用）: chatcarlo/trajectory.py
- 線量マップ換算・非物理的最大値の警告: chatcarlo/diagnostics.py
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .dose_coefficients import h_star_10_per_fluence
from .geometry import Geometry
from .materials import density, material_groups, mu_en_rho, mu_rho_parts
from .physics import (isotropic_direction, sample_compton_bound,
                       sample_fluorescence, sample_rayleigh_cos_theta,
                       sample_rayleigh_element, scatter_direction)
from .source import photon_count_through_field, sample_source_photons
from .tally import VoxelGrid, accumulate_track_length
from .trajectory import TrajectoryRecorder


def _mu_and_parts_batch(materials: np.ndarray, energies: np.ndarray):
    """線減弱係数μ[1/cm]と、光電/コンプトン/レイリー内訳の質量減弱係数[cm²/g]を、
    材料グループごとに`mu_rho_parts`を1回だけ呼んでまとめて求める。

    μ(自由行程判定用)と内訳(相互作用種別抽選用)は元々別々に
    `linear_mu`/`mu_rho_parts`を呼んでおり、同じ(材料,エネルギー)組の
    断面積テーブル参照を1ラウンドに2回行っていた（docs/plan_transport_speedup.md
    Phase 1で実測: xraylib呼び出し削減後もテーブル参照コール自体の固定
    オーバーヘッドが支配的になったため、この重複を解消した）。μは内訳の和に
    密度を掛けるだけなので、常に元の`linear_mu`と同一の値になる。
    """
    n = len(materials)
    mu = np.zeros(n)
    photo = np.zeros(n)
    compt = np.zeros(n)
    rayl = np.zeros(n)
    for name, m in material_groups(materials):
        parts = mu_rho_parts(name, energies[m])
        photo[m] = parts["photoelectric"]
        compt[m] = parts["compton"]
        rayl[m] = parts["rayleigh"]
        mu[m] = (parts["photoelectric"] + parts["compton"] + parts["rayleigh"]) * density(name)
    return mu, {"photoelectric": photo, "compton": compt, "rayleigh": rayl}


def _mu_en_linear_batch(materials: np.ndarray, energies: np.ndarray) -> np.ndarray:
    """μen = (μen/ρ)·ρ [1/cm] — カーマtrack-length estimator用の線減弱係数。"""
    mu_en = np.zeros(len(materials))
    for name, m in material_groups(materials):
        mu_en[m] = mu_en_rho(name, energies[m]) * density(name)
    return mu_en


def _deposit(energy_deposited: dict, mat_arr: np.ndarray, e_arr: np.ndarray) -> None:
    for name, m in material_groups(mat_arr):
        energy_deposited[name] = energy_deposited.get(name, 0.0) + float(np.sum(e_arr[m]))


def _interaction_probabilities_from_parts(parts_full: dict, mask: np.ndarray):
    """相互作用点ごとの(光電確率, 光電+コンプトン累積確率)。残りがレイリー。

    `_mu_and_parts_batch`が全生存光子について既に求めた内訳をmaskで
    切り出すだけで、断面積テーブルへは再アクセスしない（同一(材料,エネルギー)
    組に対する値は決定的に同じなので、フルバッチ計算からの切り出しと
    サブセットのみの独立再計算は数値的に厳密に一致する）。
    """
    photo = parts_full["photoelectric"][mask]
    compt = parts_full["compton"][mask]
    rayl = parts_full["rayleigh"][mask]
    tot = photo + compt + rayl
    tot = np.where(tot > 0, tot, 1.0)
    return photo / tot, compt / tot


@dataclass
class BatchResult:
    n_scatter: np.ndarray       # (N,) int — 相互作用回数（吸収前含む）
    absorbed: np.ndarray        # (N,) bool — 光電吸収で消滅したか（蛍光放出時はFalse）
    escaped: np.ndarray         # (N,) bool — 相互作用なしで世界境界を脱出したか
    final_energy: np.ndarray    # (N,) keV
    energy_deposited: dict = field(default_factory=dict)  # 材料名 -> keV
    n_fluorescence: int = 0     # K殻蛍光X線を放出したイベント数


def transport_photons(pos: np.ndarray, dirv: np.ndarray, energy: np.ndarray,
                       geometry: Geometry, rng: np.random.Generator,
                       grid: VoxelGrid | None = None,
                       recorder: TrajectoryRecorder | None = None,
                       tally_rng: np.random.Generator | None = None,
                       fluorescence_enabled: bool = True) -> BatchResult:
    """光源サンプリングとは独立な輸送カーネル本体（テストで直接叩ける）。

    pos/dirv/energy は呼び出し側の配列を破壊的に更新する。
    grid を渡すと、各飛行区間ごとにカーマのtrack-length estimatorを
    ボクセルグリッドへ積算する（chatcarlo/tally.py参照）。タリーの層化
    サンプリングにはtally_rng（未指定ならrngからspawnで決定的に導出）を使う。
    spawnは輸送の乱数列を消費しないため、grid有無で輸送結果（吸収/脱出・
    相互作用サンプリング）は同一seedならビット一致のまま変わらない。
    recorder を渡すと、各飛行区間を可視化用に記録する（既定Noneで無効、
    乱数を一切消費しないため同一seedでの輸送結果に影響しない）。
    fluorescence_enabled=True（既定）では光電吸収イベントでK殻蛍光X線の
    放出を抽選し（chatcarlo/physics.py の `sample_fluorescence`）、放出時は
    光子を消滅させず蛍光線エネルギー・等方方向で輸送を継続する
    （docs/plan_fluorescence.md参照）。Falseなら従来どおり全量その場で
    局所吸収する。
    """
    if grid is not None and tally_rng is None:
        tally_rng = rng.spawn(1)[0]
    n = pos.shape[0]
    alive = np.ones(n, dtype=bool)
    tau = -np.log(rng.random(n))
    n_scatter = np.zeros(n, dtype=int)
    absorbed = np.zeros(n, dtype=bool)
    escaped = np.zeros(n, dtype=bool)
    energy_deposited: dict = {}
    n_fluorescence = 0

    while np.any(alive):
        idx = np.where(alive)[0]
        o, d, e = pos[idx], dirv[idx], energy[idx]
        mat = geometry.material_at(o)
        mu, parts_full = _mu_and_parts_batch(mat, e)
        t_boundary, escape = geometry.next_boundary(o, d)
        mu_safe = np.where(mu > 0, mu, 1e-30)
        tau_to_boundary = mu * t_boundary
        will_interact = tau[idx] < tau_to_boundary

        ds = np.where(will_interact, tau[idx] / mu_safe, t_boundary)
        ends = o + d * ds[:, None]

        if grid is not None:
            mu_en_linear = _mu_en_linear_batch(mat, e)
            accumulate_track_length(grid.kerma_keV, grid, o, d, ds, e * mu_en_linear, tally_rng)
            accumulate_track_length(grid.h10_track_pSv_cm3, grid, o, d, ds,
                                     h_star_10_per_fluence(e), tally_rng)

        pos[idx] = ends

        noninteract = ~will_interact
        gidx = idx[noninteract]
        tau[gidx] -= tau_to_boundary[noninteract]
        pos[gidx] += dirv[gidx] * 1e-6
        esc_now = idx[noninteract & escape]
        alive[esc_now] = False
        escaped[esc_now] = True

        interact = will_interact
        iidx = idx[interact]
        if len(iidx) > 0:
            mat_i = mat[interact]
            e_i = e[interact]
            r_type = rng.random(len(iidx))
            p_photo, p_compt = _interaction_probabilities_from_parts(parts_full, interact)

            is_photo = r_type < p_photo
            is_compt = (~is_photo) & (r_type < p_photo + p_compt)
            is_rayl = (~is_photo) & (~is_compt)

            photo_idx = iidx[is_photo]
            is_fluor = np.zeros(len(iidx), dtype=bool)  # is_photoと同じ形（recorder用）
            if len(photo_idx) > 0:
                mat_p = mat_i[is_photo]
                e_p = e_i[is_photo]
                if fluorescence_enabled:
                    emit, e_line = sample_fluorescence(mat_p, e_p, rng)
                else:
                    emit = np.zeros(len(photo_idx), dtype=bool)
                    e_line = np.zeros(len(photo_idx))

                _deposit(energy_deposited, mat_p, np.where(emit, e_p - e_line, e_p))
                n_scatter[photo_idx] += 1

                no_emit = photo_idx[~emit]
                alive[no_emit] = False
                absorbed[no_emit] = True

                emit_idx = photo_idx[emit]
                if len(emit_idx) > 0:
                    n_fluorescence += len(emit_idx)
                    energy[emit_idx] = e_line[emit]
                    dirv[emit_idx] = isotropic_direction(len(emit_idx), rng)
                    tau[emit_idx] = -np.log(rng.random(len(emit_idx)))

                photo_positions = np.where(is_photo)[0]
                is_fluor[photo_positions] = emit

            compt_idx = iidx[is_compt]
            if len(compt_idx) > 0:
                e_c = e_i[is_compt]
                eps, cos_theta = sample_compton_bound(mat_i[is_compt], e_c, rng)
                e_new = e_c * eps
                _deposit(energy_deposited, mat_i[is_compt], e_c - e_new)
                dirv[compt_idx] = scatter_direction(dirv[compt_idx], cos_theta, rng)
                energy[compt_idx] = e_new
                tau[compt_idx] = -np.log(rng.random(len(compt_idx)))
                n_scatter[compt_idx] += 1

            rayl_idx = iidx[is_rayl]
            if len(rayl_idx) > 0:
                z_r = sample_rayleigh_element(mat_i[is_rayl], e_i[is_rayl], rng)
                cos_theta = sample_rayleigh_cos_theta(z_r, e_i[is_rayl], rng)
                dirv[rayl_idx] = scatter_direction(dirv[rayl_idx], cos_theta, rng)
                tau[rayl_idx] = -np.log(rng.random(len(rayl_idx)))
                n_scatter[rayl_idx] += 1

        if recorder is not None:
            event = np.full(len(idx), "boundary", dtype=object)
            event[noninteract & escape] = "escape"
            if len(iidx) > 0:
                interact_positions = np.where(interact)[0]
                photo_positions_full = interact_positions[is_photo]
                event[photo_positions_full] = "photoelectric"
                event[photo_positions_full[is_fluor[is_photo]]] = "fluorescence"
                event[interact_positions[is_compt]] = "compton"
                event[interact_positions[is_rayl]] = "rayleigh"
            recorder.record(idx, o, ends, e, event)

    return BatchResult(n_scatter=n_scatter, absorbed=absorbed, escaped=escaped,
                        final_energy=energy, energy_deposited=energy_deposited,
                        n_fluorescence=n_fluorescence)


@dataclass
class TransportResult:
    n_histories: int
    energy_deposited_MeV: dict
    fraction_absorbed: float
    fraction_escaped: float
    mean_scatter_events: float
    grid: VoxelGrid | None = None
    # 絶対値換算係数（per-history値×これ=実線量）。mas指定時は照射野を通過する
    # 実光子数、ctdi_vol_mGy指定時はCTDIファントム校正による実効光子数。
    n_photons_real: float | None = None
    n_fluorescence: int = 0


def run_transport(scene, n_histories: int = 100_000, seed: int | None = None,
                   batch_size: int = 200_000, dose_grid: bool = False,
                   grid_resolution_cm: float = 5.0) -> TransportResult:
    rng = np.random.default_rng(seed)
    src = scene.raw["source"]
    geometry = Geometry(scene.raw["geometry"])
    grid = VoxelGrid.from_bbox(geometry.bbox_min, geometry.bbox_max, grid_resolution_cm) if dose_grid else None
    fluorescence_enabled = scene.raw.get("physics", {}).get("fluorescence", True)

    energy_deposited: dict = {}
    n_absorbed = 0
    n_escaped = 0
    scatter_sum = 0
    n_fluorescence = 0
    remaining = n_histories
    while remaining > 0:
        n = min(batch_size, remaining)
        remaining -= n
        pos, dirv, energy = sample_source_photons(src, n, rng)
        result = transport_photons(pos, dirv, energy, geometry, rng, grid=grid,
                                    fluorescence_enabled=fluorescence_enabled)
        for name, e_keV in result.energy_deposited.items():
            energy_deposited[name] = energy_deposited.get(name, 0.0) + e_keV
        n_fluorescence += result.n_fluorescence
        n_absorbed += int(np.sum(result.absorbed))
        n_escaped += int(np.sum(result.escaped))
        scatter_sum += int(np.sum(result.n_scatter))

    # 絶対線量校正: CTDIvol基準（CT向け、実測に装置特性が折り込まれ汎用性が高い）
    # が指定されていればそちらを優先。なければmAs+SpekPyフルエンス基準。
    if src.get("ctdi_vol_mGy") is not None:
        from .ctdi import effective_histories_from_ctdi
        n_photons_real = effective_histories_from_ctdi(src, seed=seed)
    elif src.get("mas") is not None:
        n_photons_real = photon_count_through_field(src)
    else:
        n_photons_real = None

    return TransportResult(
        n_histories=n_histories,
        energy_deposited_MeV={k: v / 1000.0 for k, v in energy_deposited.items()},
        fraction_absorbed=n_absorbed / n_histories,
        fraction_escaped=n_escaped / n_histories,
        mean_scatter_events=scatter_sum / n_histories,
        grid=grid,
        n_photons_real=n_photons_real,
        n_fluorescence=n_fluorescence,
    )
