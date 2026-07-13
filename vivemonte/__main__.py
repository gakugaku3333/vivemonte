"""CLI — AI（Claude Code）が非対話で叩ける入口。

  python -m vivemonte validate <scene.yaml>
  python -m vivemonte preview  <scene.yaml> [-o out.html]
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
