"""CLI — AI（Claude Code）が非対話で叩ける入口。

  python -m vivemonte validate <scene.yaml>
  python -m vivemonte preview  <scene.yaml> [-o out.html]
  python -m vivemonte trace    <scene.yaml> [-n 200] [--seed 42] [-o out.html]
  python -m vivemonte xs <材料...> [--emin 10] [--emax 150] [-o out.png]
"""
from __future__ import annotations

import argparse
import sys


def cmd_validate(args) -> int:
    from .scene import load_scene
    scene = load_scene(args.scene)
    for w in scene.warnings:
        print(f"[警告] {w}")
    if scene.ok:
        print(f"OK: {args.scene} は有効なシーンです（物体 {len(scene.raw['geometry'])} 個）")
        return 0
    for e in scene.errors:
        print(f"[エラー] {e}", file=sys.stderr)
    return 1


def cmd_preview(args) -> int:
    from .preview import write_html
    from .scene import load_scene
    scene = load_scene(args.scene)
    if not scene.ok:
        for e in scene.errors:
            print(f"[エラー] {e}", file=sys.stderr)
        return 1
    for w in scene.warnings:
        print(f"[警告] {w}")
    out = args.out or args.scene.rsplit(".", 1)[0] + "_preview.html"
    write_html(scene, out, title=args.title or f"viveMonte — {args.scene}")
    print(f"プレビューを書き出しました: {out}")
    return 0


_TRACE_MAX_N = 2000


def cmd_trace(args) -> int:
    import numpy as np

    from .geometry import Geometry
    from .preview import write_html
    from .scene import load_scene
    from .transport import (TrajectoryRecorder, sample_source_photons,
                             transport_photons, trajectories_to_json)

    scene = load_scene(args.scene)
    if not scene.ok:
        for e in scene.errors:
            print(f"[エラー] {e}", file=sys.stderr)
        return 1
    for w in scene.warnings:
        print(f"[警告] {w}")

    n = int(args.n)
    if n > _TRACE_MAX_N:
        print(f"[警告] -n={n} は大きすぎるため{_TRACE_MAX_N}にクランプします"
              "（HTML描画性能のため）")
        n = _TRACE_MAX_N

    rng = np.random.default_rng(args.seed)
    geometry = Geometry(scene.raw["geometry"])
    pos, dirv, energy = sample_source_photons(scene.raw["source"], n, rng)
    recorder = TrajectoryRecorder()
    transport_photons(pos, dirv, energy, geometry, rng, recorder=recorder)
    trajectories = trajectories_to_json(recorder)

    out = args.out or args.scene.rsplit(".", 1)[0] + "_trace.html"
    write_html(scene, out, title=args.title or f"viveMonte trace — {args.scene}",
               trajectories=trajectories)
    print(f"光子{n}個の軌跡を書き出しました: {out}")
    return 0


def cmd_run(args) -> int:
    from .scene import load_scene
    from .transport import dose_map_Gy, run_transport
    from .geometry import Geometry

    scene = load_scene(args.scene)
    if not scene.ok:
        for e in scene.errors:
            print(f"[エラー] {e}", file=sys.stderr)
        return 1
    for w in scene.warnings:
        print(f"[警告] {w}")

    result = run_transport(scene, n_histories=int(args.n_histories), seed=args.seed,
                            dose_grid=args.dose_grid, grid_resolution_cm=args.resolution)
    print(f"histories: {result.n_histories:,}")
    print(f"吸収（光電）割合: {result.fraction_absorbed:.4f}")
    print(f"脱出割合: {result.fraction_escaped:.4f}")
    print(f"平均相互作用回数/光子: {result.mean_scatter_events:.4f}")
    print("材料別吸収エネルギー [MeV/history合計]:")
    for name, e_mev in sorted(result.energy_deposited_MeV.items(), key=lambda kv: -kv[1]):
        print(f"  {name}: {e_mev:.6g}")

    if result.n_photons_real is not None:
        print(f"光子数校正: mAs={scene.raw['source']['mas']:g} で照射野を通過する実光子数 "
              f"= {result.n_photons_real:.6g}")

    if args.dose_grid:
        import numpy as np
        geometry = Geometry(scene.raw["geometry"])
        dose = dose_map_Gy(result.grid, geometry)
        h10 = result.grid.h10_map_pSv()
        n_histories = result.n_histories
        # scene.source.mas未指定時は1historyあたりの値のみ（相対値）を出力する。
        # mas指定時は実光子数でスケールした絶対値[Gy]・[pSv]も出す。
        # dose_per_history/h10_per_historyは既にn_historiesで割って「1光子あたり」に
        # なっているので、ここでの係数は実光子数そのもの（n_historiesでは割らない）。
        scale = result.n_photons_real if result.n_photons_real is not None else None
        dose_per_history = dose / n_histories
        h10_per_history = h10 / n_histories
        print(f"線量グリッド: shape={result.grid.shape}, "
              f"resolution={result.grid.voxel_size_cm}cm")
        print(f"  最大吸収線量 [Gy/history]: {dose_per_history.max():.6g}")
        print(f"  総カーマ: {result.grid.total_kerma_MeV():.6g} MeV "
              f"({result.grid.total_kerma_MeV() / n_histories * 1000:.6g} keV/history)")
        print(f"  最大H*(10) [pSv/history]: {h10_per_history.max():.6g}")
        if scale is not None:
            print(f"  最大吸収線量 [Gy]（mAs校正済み）: {(dose_per_history.max() * scale):.6g}")
            print(f"  最大H*(10) [pSv]（mAs校正済み）: {(h10_per_history.max() * scale):.6g}")
        if args.dose_out:
            save_kwargs = dict(
                dose_per_history_Gy=dose_per_history,
                h10_per_history_pSv=h10_per_history,
                kerma_keV=result.grid.kerma_keV,
                h10_track_pSv_cm3=result.grid.h10_track_pSv_cm3,
                origin_cm=result.grid.origin_cm,
                voxel_size_cm=result.grid.voxel_size_cm,
                shape=np.array(result.grid.shape),
            )
            if scale is not None:
                save_kwargs["dose_Gy"] = dose_per_history * scale
                save_kwargs["h10_pSv"] = h10_per_history * scale
            np.savez(args.dose_out, **save_kwargs)
            print(f"線量グリッドを書き出しました: {args.dose_out}")
    return 0


def cmd_xs(args) -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    from .materials import mu_en_rho, mu_rho, resolve

    e = np.logspace(np.log10(args.emin), np.log10(args.emax), 200)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    for m in args.materials:
        name, rho, _ = resolve(m)
        ax1.loglog(e, mu_rho(m, e), label=f"{m} (ρ={rho:.3g})")
        ax2.loglog(e, mu_en_rho(m, e), label=m)
    ax1.set_title("Mass attenuation μ/ρ  [cm²/g]")
    ax2.set_title("Mass energy-absorption μen/ρ  [cm²/g]")
    for ax in (ax1, ax2):
        ax.set_xlabel("Photon energy [keV]")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("viveMonte cross-section data (xraylib / EPDL, NIST-consistent)")
    fig.tight_layout()
    out = args.out or "xs.png"
    fig.savefig(out, dpi=140)
    print(f"断面積プロットを書き出しました: {out}")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(prog="vivemonte")
    sub = p.add_subparsers(dest="cmd", required=True)

    pv = sub.add_parser("validate", help="scene.yaml を検証")
    pv.add_argument("scene")
    pv.set_defaults(func=cmd_validate)

    pp = sub.add_parser("preview", help="ジオメトリープレビューHTMLを生成")
    pp.add_argument("scene")
    pp.add_argument("-o", "--out")
    pp.add_argument("--title")
    pp.set_defaults(func=cmd_preview)

    pt = sub.add_parser("trace", help="光子軌跡を3D可視化するHTMLを生成（小history）")
    pt.add_argument("scene")
    pt.add_argument("-n", type=int, default=200, help="軌跡を記録する光子数（既定200、上限2000）")
    pt.add_argument("--seed", type=int, default=None)
    pt.add_argument("-o", "--out")
    pt.add_argument("--title")
    pt.set_defaults(func=cmd_trace)

    pr = sub.add_parser("run", help="光子輸送を実行")
    pr.add_argument("scene")
    pr.add_argument("-n", "--n-histories", type=float, default=1e5)
    pr.add_argument("--seed", type=int, default=None)
    pr.add_argument("--dose-grid", action="store_true", help="ボクセル吸収線量タリーを有効化")
    pr.add_argument("--resolution", type=float, default=5.0, help="線量グリッド解像度[cm]（既定5cm）")
    pr.add_argument("--dose-out", help="線量グリッドを.npzに書き出すパス")
    pr.set_defaults(func=cmd_run)

    px = sub.add_parser("xs", help="断面積カーブを描画")
    px.add_argument("materials", nargs="+")
    px.add_argument("--emin", type=float, default=10.0)
    px.add_argument("--emax", type=float, default=150.0)
    px.add_argument("-o", "--out")
    px.set_defaults(func=cmd_xs)

    args = p.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
