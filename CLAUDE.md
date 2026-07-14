# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 報告時のルール

作業完了時にユーザーへ報告する際は、進捗の説明だけで終わらせず、成果物そのものを
提示する（生成したHTML/PNG/npzファイルを`SendUserFile`等で送る、重要な数値結果は
本文に具体的に書く、など）。「〜を実装しました」で終わらせて中身を見せないのは
避けること。

## What this is

viveMonte is a Monte Carlo photon transport code for diagnostic X-ray energies (10–150 keV). Scenes are declared
in `scene.yaml` — the system is designed for an AI (Claude Code) to write and iterate on scene files declaratively,
then drive validation/preview/run non-interactively via the CLI. Current status is research/education only: doses
and H*(10) must be cross-checked against an established code (PHITS) before being used for real patient-dose or
shielding decisions (see the warning banner in [README.md](README.md)).

## Commands

```bash
# setup (venv is project-local, per the parent Projects/CLAUDE.md rule — don't pip install globally)
python3 -m venv .venv
.venv/bin/pip install numpy pyyaml matplotlib xraylib pytest spekpy

# validate a scene (physical sanity checks, not just schema)
.venv/bin/python -m vivemonte validate examples/chest_room.yaml

# 3D geometry preview -> self-contained HTML (no external deps)
.venv/bin/python -m vivemonte preview examples/chest_room.yaml -o preview.html

# cross-section curves
.venv/bin/python -m vivemonte xs water bone lead -o xs.png

# run transport (prints per-material absorbed energy, absorbed/escaped fractions)
.venv/bin/python -m vivemonte run examples/chest_room.yaml -n 1e6 --seed 42

# same, plus voxel absorbed-dose/H*(10) tally written to .npz
.venv/bin/python -m vivemonte run examples/chest_room.yaml -n 1e6 --seed 42 \
    --dose-grid --resolution 5 --dose-out dose.npz

# photon trajectory 3D visualization (small n; overlays onto the preview HTML template)
.venv/bin/python -m vivemonte trace examples/chest_room.yaml -n 200 --seed 42 -o trace.html

# cross-section slices through a dose/H*(10) map (default: 3 planes through the max-value voxel)
.venv/bin/python -m vivemonte plot dose.npz --scene examples/chest_room.yaml -o maps.png

# tests (spot-checked against published NIST reference values)
.venv/bin/python -m pytest tests/ -q
```

Run a single test: `.venv/bin/python -m pytest tests/test_transport.py::test_name -q`

There's also `.claude/skills/vive-check/`, a gated workflow skill that runs the four CLI steps in order
(geometry preview → trajectory preview → full run → results) with human approval at each gate. Invoke it for
"walk through the scene with me" style requests rather than chaining the raw CLI calls yourself.

When the user asks for a simulation but no scene.yaml exists yet (or the request is vague), do NOT start writing
a scene directly — invoke `.claude/skills/vive-interview/` first. It elicits the requirements in stages
(purpose → exposure parameters → geometry → run settings) via AskUserQuestion, confirming intent and pinning down
ambiguities before drafting scene.yaml, then hands off to vive-check.

## Architecture

**Transport is not voxelized.** Geometry stays as analytic primitives (box/cylinder/sphere in
[geometry.py](vivemonte/geometry.py)); each photon steps to the next material boundary by computing an analytic
ray/primitive intersection distance ([transport.py](vivemonte/transport.py)). This is "analytic surface tracking,"
not Woodcock delta-tracking — since each segment has one homogeneous material, μ is constant along a step, so no
virtual-collision rejection is needed. This scales well for room-size scenes with thin shielding without a
voxel-resolution/memory tradeoff. Overlapping objects resolve by list order (later wins); open space not inside
any object is `background` (default air).

**Module layout around the kernel**: [transport.py](vivemonte/transport.py) is only the transport loop +
`run_transport`. Spectrum generation (SpekPy/Kramers, heel off-axis spectra) lives in
[spectrum.py](vivemonte/spectrum.py); source/field sampling and the mAs photon-count calibration in
[source.py](vivemonte/source.py); interaction angle/energy sampling in [physics.py](vivemonte/physics.py);
trajectory recording for `trace` in [trajectory.py](vivemonte/trajectory.py); dose-map conversion and the
non-physical-max warnings in [diagnostics.py](vivemonte/diagnostics.py).

**Physics**: photoelectric / Compton (Klein-Nishina with Kahn rejection sampling) / Rayleigh (atomic form factor
F(Z,q) via `xraylib.FF_Rayl`, compounds sampled by mass-fraction-weighted element pick before the angular
distribution). Electron range is neglected (kerma approximation — local absorption at the interaction point).
[tests/test_transport.py](tests/test_transport.py) checks primary transmission against the analytic Beer-Lambert
law (`exp(-μt)`).

**Dose/H*(10) tallying is a separate concern from transport.** [tally.py](vivemonte/tally.py)'s `VoxelGrid` lays a
uniform grid independently of the transport geometry, purely for scoring. Two independent estimators are
cross-validated against each other in [tests/test_tally.py](tests/test_tally.py): a collision estimator
(`energy_deposited`, scored at interaction points) and a track-length kerma estimator (path-integral over the
grid). H*(10) is a fluence-based protection quantity (different from kerma), computed by normalizing the
track-length integral by voxel volume (`VoxelGrid.h10_map_pSv`).

**Units and calibration**: relative output is `Gy/history` / `pSv/history`. When `scene.yaml`'s `source.mas` is
set, `photon_count_through_field` (in source.py) uses SpekPy's absolute fluence to get the real photon count
through the field, and per-history values are scaled by that count (not divided again by `n_histories` — see the
mAs double-division bug writeup in [docs/lessons_learned.md](docs/lessons_learned.md) if touching this path).

**Cross-section/dose-coefficient data provenance** (do not "improve" these from memory — see lessons learned):
| quantity | source | used for |
|---|---|---|
| μ/ρ, photoelectric/Compton/Rayleigh split | `xraylib` (EPDL-based, matches NIST XCOM) | transport free-path/interaction sampling |
| μen/ρ (mass energy-absorption coefficient) | NIST XAAMDI, bundled CSV (`vivemonte/data/nist_xaamdi/`) | kerma/absorbed-dose tally |
| h*(10)/Φ (ambient dose equivalent conversion) | ICRP Publication 74 / ICRU Report 57, bundled CSV (`vivemonte/data/h_star_10/`) | H*(10) tally |

`xraylib`'s `CS_Energy` diverges from NIST-published μen/ρ by up to ~17% and must not be used for dose — this is
covered by a regression test in `tests/test_materials.py`. Re-fetch data with `scripts/fetch_nist_xaamdi.py` /
`scripts/fetch_h_star_10.py`, never by hand-typing values.

## Known sharp edges (read before trusting "max dose"/"max H*(10)" output)

`vivemonte run --dose-grid`'s reported max absorbed dose / max H*(10) frequently lands in a non-physical spot: a
background (air) voxel near the 1/r² point-source singularity, or an air voxel just outside a material boundary
due to backscatter (this gets worse, not better, with finer grid resolution — it's a boundary effect, not
discretization error). The CLI detects and prints a `[警告]` when this happens
(`background_medium_warning`/`near_source_air_warning` in diagnostics.py). For a real exposure-point estimate (patient
surface, operator position, etc.), lay a fine grid directly at that position rather than trusting the global max.
Full writeup: [docs/lessons_learned.md](docs/lessons_learned.md).

## Scene files

`scene.yaml`: cm units, z-axis up, floor at z=0. See [vivemonte/scene.py](vivemonte/scene.py) for the validator
(`load_scene`/`validate_scene`) — it's designed to produce actionable errors (`geometry[2].size_cm: ...`) for an AI
self-correction loop, plus non-fatal physics-sanity warnings (e.g. filtration below legal minimum, source position
inside a solid). [examples/chest_room.yaml](examples/chest_room.yaml) is the canonical worked example (standing
chest X-ray room).
