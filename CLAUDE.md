# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 報告時のルール

作業完了時にユーザーへ報告する際は、進捗の説明だけで終わらせず、成果物そのものを
提示する（生成したHTML/PNG/npzファイルを`SendUserFile`等で送る、重要な数値結果は
本文に具体的に書く、など）。「〜を実装しました」で終わらせて中身を見せないのは
避けること。

## What this is

ChatCarlo is a Monte Carlo photon transport code for diagnostic X-ray energies (10–150 keV). Scenes are declared
in `scene.yaml` — the system is designed for an AI (Claude Code) to write and iterate on scene files declaratively,
then drive validation/preview/run non-interactively via the CLI. Current status is research/education only: doses
and H*(10) must be cross-checked against an established code (EGS5 — PHITS was considered but was not adopted;
see [docs/plan_egs5_crosscheck.md](docs/plan_egs5_crosscheck.md)) before being used for real patient-dose or
shielding decisions (see the warning banner in [README.md](README.md)).

## Commands

```bash
# setup (venv is project-local, per the parent Projects/CLAUDE.md rule — don't pip install globally)
python3 -m venv .venv
.venv/bin/pip install numpy pyyaml matplotlib xraylib pytest spekpy scipy

# validate a scene (physical sanity checks, not just schema)
.venv/bin/python -m chatcarlo validate examples/chest_room.yaml

# 3D geometry preview -> self-contained HTML (no external deps)
.venv/bin/python -m chatcarlo preview examples/chest_room.yaml -o preview.html

# cross-section curves
.venv/bin/python -m chatcarlo xs water bone lead -o xs.png

# run transport (prints per-material absorbed energy, absorbed/escaped fractions)
.venv/bin/python -m chatcarlo run examples/chest_room.yaml -n 1e6 --seed 42

# same, plus voxel absorbed-dose/H*(10) tally written to .npz
.venv/bin/python -m chatcarlo run examples/chest_room.yaml -n 1e6 --seed 42 \
    --dose-grid --resolution 5 --dose-out dose.npz

# photon trajectory 3D visualization (small n; overlays onto the preview HTML template)
.venv/bin/python -m chatcarlo trace examples/chest_room.yaml -n 200 --seed 42 -o trace.html

# cross-section slices through a dose/H*(10) map (default: 3 planes through the max-value voxel)
.venv/bin/python -m chatcarlo plot dose.npz --scene examples/chest_room.yaml -o maps.png

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
[geometry.py](chatcarlo/geometry.py)); each photon steps to the next material boundary by computing an analytic
ray/primitive intersection distance ([transport.py](chatcarlo/transport.py)). This is "analytic surface tracking,"
not Woodcock delta-tracking — since each segment has one homogeneous material, μ is constant along a step, so no
virtual-collision rejection is needed. This scales well for room-size scenes with thin shielding without a
voxel-resolution/memory tradeoff. Overlapping objects resolve by list order (later wins); open space not inside
any object is `background` (default air).

**Module layout around the kernel**: [transport.py](chatcarlo/transport.py) is only the transport loop +
`run_transport`. Spectrum generation (SpekPy/Kramers, heel off-axis spectra) lives in
[spectrum.py](chatcarlo/spectrum.py); source/field sampling and the mAs photon-count calibration in
[source.py](chatcarlo/source.py); interaction angle/energy sampling in [physics.py](chatcarlo/physics.py);
trajectory recording for `trace` in [trajectory.py](chatcarlo/trajectory.py); dose-map conversion and the
non-physical-max warnings in [diagnostics.py](chatcarlo/diagnostics.py).

**Physics**: photoelectric / Compton (bound Compton — free-electron Klein-Nishina via Kahn rejection sampling,
then an additional S(Z,q)/Z rejection from the incoherent scattering function via `xraylib.SF_Compt`; compounds
sampled by mass-fraction-weighted element pick, same pattern as Rayleigh, before the angular distribution) /
Rayleigh (atomic form factor F(Z,q) via `xraylib.FF_Rayl`, compounds sampled by mass-fraction-weighted element
pick before the angular distribution). Electron range is neglected (kerma approximation — local absorption at
the interaction point), except that photoelectric absorption samples K-shell fluorescence emission
(`sample_fluorescence` in [physics.py](chatcarlo/physics.py); K-shell only, no cascade/L-shell, line energies
below 5 keV are absorbed locally instead of emitted) — when emitted, the photon continues transport at the
fluorescence line energy with an isotropic direction rather than being annihilated. Controlled by
`physics.fluorescence` in scene.yaml (default `true`); toggling it off reproduces the pre-fluorescence local-absorption
behavior. See [docs/plan_fluorescence.md](docs/plan_fluorescence.md) for the design rationale and verification.
[tests/test_transport.py](tests/test_transport.py) checks primary transmission against the analytic Beer-Lambert
law (`exp(-μt)`); [tests/test_fluorescence.py](tests/test_fluorescence.py) checks K-edge data against xraylib,
energy conservation with fluorescence on/off, and the emission rate against the analytic K-shell-fraction×ω_K
expectation.

**Dose/H*(10) tallying is a separate concern from transport.** [tally.py](chatcarlo/tally.py)'s `VoxelGrid` lays a
uniform grid independently of the transport geometry, purely for scoring. Two independent estimators are
cross-validated against each other in [tests/test_tally.py](tests/test_tally.py): a collision estimator
(`energy_deposited`, scored at interaction points) and a track-length kerma estimator (path-integral over the
grid). H*(10) is a fluence-based protection quantity (different from kerma), computed by normalizing the
track-length integral by voxel volume (`VoxelGrid.h10_map_pSv`). `accumulate_track_length` splits each flight
segment into substeps and scores a **stratified random point within each substep** (not the substep midpoint) —
this makes the spatial-binning step an unbiased estimator regardless of substep length, which matters when many
segments start exactly on a voxel boundary (e.g. a `field.shape: parallel` beam entering a phantom face); see
lessons_learned for the bug this replaced. Before K-shell fluorescence was modeled, the two estimators disagreed
by design in high-Z materials — the collision estimator deposited the full photoelectric energy locally while
NIST μen/ρ (used by the track-length estimator) already subtracts the mean fluorescence escape fraction. Modeling
fluorescence brought the two into much closer agreement for lead (spot-checked with a 100 keV beam into a thick
lead slab: track-length/collision ratio improved from ~0.38 without fluorescence to ~0.92 with it — see
[docs/plan_fluorescence.md](docs/plan_fluorescence.md) for the verification script and numbers).

**Units and calibration**: relative output is `Gy/history` / `pSv/history`. When `scene.yaml`'s `source.mas` is
set, `photon_count_through_field` (in source.py) uses SpekPy's absolute fluence to get the real photon count
through the field, and per-history values are scaled by that count (not divided again by `n_histories` — see the
mAs double-division bug writeup in [docs/lessons_learned.md](docs/lessons_learned.md) if touching this path).

**Cross-section/dose-coefficient data provenance** (do not "improve" these from memory — see lessons learned):
| quantity | source | used for |
|---|---|---|
| μ/ρ, photoelectric/Compton/Rayleigh split | `xraylib` (EPDL-based, matches NIST XCOM) | transport free-path/interaction sampling |
| μen/ρ (mass energy-absorption coefficient) | NIST XAAMDI, bundled CSV (`chatcarlo/data/nist_xaamdi/`) | kerma/absorbed-dose tally |
| h*(10)/Φ (ambient dose equivalent conversion) | ICRP Publication 74 / ICRU Report 57, bundled CSV (`chatcarlo/data/h_star_10/`) | H*(10) tally |

`xraylib`'s `CS_Energy` diverges from NIST-published μen/ρ by up to ~17% and must not be used for dose — this is
covered by a regression test in `tests/test_materials.py`. Re-fetch data with `scripts/fetch_nist_xaamdi.py` /
`scripts/fetch_h_star_10.py`, never by hand-typing values.

## Known sharp edges (read before trusting "max dose"/"max H*(10)" output)

`chatcarlo run --dose-grid`'s reported max absorbed dose / max H*(10) can still land in a non-physical spot: a
background (air) voxel near the 1/r² point-source singularity, or an air voxel just outside a material boundary
due to backscatter. The CLI detects and prints a `[警告]` when this happens
(`background_medium_warning`/`near_source_air_warning` in diagnostics.py). For a real exposure-point estimate (patient
surface, operator position, etc.), lay a fine grid directly at that position rather than trusting the global max.

A related but distinct effect: at fine resolution (≲0.5cm) the reported max can still grow as resolution gets
finer. This is **not** the same bug as a since-fixed systematic bias where flight segments starting exactly on a
voxel boundary (all photons of a `field.shape: parallel` beam, for instance) were under-scored in that boundary
voxel — that one was a real bug in the track-length tally's substep scoring and is fixed (see
`accumulate_track_length` above). The residual fine-resolution growth is believed to be extreme-value statistics
from `max_substeps` clamping on long segments, not yet independently re-verified below ~1cm. Both warnings above
only fire when the max lands on `background`/air, so they don't catch this class of issue when both neighboring
voxels are legitimately the declared material — treat single-voxel maxima at declared-material boundaries with
the same caution. Full writeup: [docs/lessons_learned.md](docs/lessons_learned.md).

## Scene files

`scene.yaml`: cm units, z-axis up, floor at z=0. See [chatcarlo/scene.py](chatcarlo/scene.py) for the validator
(`load_scene`/`validate_scene`) — it's designed to produce actionable errors (`geometry[2].size_cm: ...`) for an AI
self-correction loop, plus non-fatal physics-sanity warnings (e.g. filtration below legal minimum, source position
inside a solid). [examples/chest_room.yaml](examples/chest_room.yaml) is the canonical worked example (standing
chest X-ray room).

**Non-clinical/research sources**: `source.kvp` (polychromatic SpekPy/Kramers spectrum) is the default and required
unless `source.spectrum` is given instead — an explicit `[{energy_keV, weight}, ...]` list (e.g. a single entry for
a monoenergetic beam), used for physics cross-checks against reference codes rather than clinical scenes. `kvp` and
`spectrum` are mutually exclusive (validated), as are `spectrum` with `mas`/`ctdi_vol_mGy`/`heel_effect` (those are
fixed to the kvp-based SpekPy path and would silently disagree with an overridden spectrum) — such scenes only
produce relative `Gy/history` output. Similarly `source.field.shape: parallel` (non-divergent beam, `size_cm` only,
no `sid_cm`) is available alongside `rect`/`cone` for the same reference-code-matching use case, with the same
`mas`/`heel_effect` restriction. See [examples/water_phantom_pdd_ocr.yaml](examples/water_phantom_pdd_ocr.yaml).
