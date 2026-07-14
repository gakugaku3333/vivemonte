"""光子軌跡レコーダー（TrajectoryRecorder / trajectories_to_json）のテスト。"""
from __future__ import annotations

import numpy as np

from vivemonte.geometry import Geometry
from vivemonte.preview import write_html
from vivemonte.scene import validate_scene
from vivemonte.source import sample_source_photons
from vivemonte.trajectory import TrajectoryRecorder, trajectories_to_json
from vivemonte.transport import transport_photons

_BASE_SOURCE = {
    "kvp": 80, "position": [0, -50, 0], "direction": [0, 1, 0],
    "field": {"size_cm": [10, 10], "sid_cm": 50}, "filtration_mm_al": 2.5,
}
_BASE_GEOMETRY = [{
    "name": "slab", "shape": "box", "material": "water",
    "center": [0, 0, 0], "size_cm": [10, 10, 10],
}]


def _slab_geometry():
    return Geometry([{
        "name": "slab", "shape": "box", "material": "water",
        "center": [0.0, 0.0, 0.0],
        "size_cm": [10.0, 100.0, 100.0],
    }])


def _make_beam(n, energy_keV):
    pos = np.tile(np.array([-30.0, 0.0, 0.0]), (n, 1))
    dirv = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    energy = np.full(n, energy_keV)
    return pos, dirv, energy


def _run(n, seed, energy_keV=60.0, recorder=None):
    geom = _slab_geometry()
    rng = np.random.default_rng(seed)
    pos, dirv, energy = _make_beam(n, energy_keV)
    result = transport_photons(pos, dirv, energy, geom, rng, recorder=recorder)
    return result


def _photon_segments(recorder: TrajectoryRecorder, photon_id: int):
    """指定photon_idの(start, end)区間を反復順（=時系列順）で取り出す。"""
    segs = []
    for ids, starts, ends in zip(recorder.photon_ids, recorder.starts, recorder.ends):
        m = ids == photon_id
        if np.any(m):
            segs.append((starts[m][0], ends[m][0]))
    return segs


def test_trajectory_continuity_single_photon():
    """区間の連続性: 各区間の始点は前区間の終点と一致する（境界ナッジ分のみ許容）。"""
    rec = TrajectoryRecorder()
    _run(n=1, seed=1, recorder=rec)
    segs = _photon_segments(rec, 0)
    assert len(segs) >= 2
    for (_, end_prev), (start_next, _) in zip(segs, segs[1:]):
        assert np.allclose(start_next, end_prev, atol=1e-5)


def test_energy_monotonic_and_conserved_except_compton():
    """エネルギー単調性: energiesは非増加。コンプトンでのみ減り、境界・レイリーでは不変。"""
    rec = TrajectoryRecorder()
    _run(n=500, seed=2, energy_keV=80.0, recorder=rec)
    trajs = trajectories_to_json(rec)
    assert trajs

    saw_compton_drop = False
    saw_unchanged = False
    for t in trajs:
        energies, events = t["energies"], t["events"]
        for i in range(len(energies) - 1):
            assert energies[i + 1] <= energies[i] + 1e-9
            if events[i] == "compton":
                assert energies[i + 1] < energies[i]
                saw_compton_drop = True
            elif events[i] in ("boundary", "rayleigh"):
                assert abs(energies[i + 1] - energies[i]) < 1e-9
                saw_unchanged = True
    assert saw_compton_drop
    assert saw_unchanged


def test_last_event_is_terminal_and_photoelectric_not_midlist():
    """イベント整合性: 最後のeventは必ずphotoelectricかescape。途中には現れない。"""
    rec = TrajectoryRecorder()
    _run(n=500, seed=3, energy_keV=60.0, recorder=rec)
    trajs = trajectories_to_json(rec)
    assert trajs
    for t in trajs:
        events = t["events"]
        assert events[-1] in ("photoelectric", "escape")
        assert "photoelectric" not in events[:-1]


def test_recorder_none_does_not_change_transport_result():
    """recorder=None の無影響: 同一seedでrecorder有無のBatchResultが完全一致する。"""
    n, seed, energy_keV = 300, 7, 60.0

    geom_a = _slab_geometry()
    rng_a = np.random.default_rng(seed)
    pos_a, dirv_a, energy_a = _make_beam(n, energy_keV)
    result_plain = transport_photons(pos_a, dirv_a, energy_a, geom_a, rng_a)

    geom_b = _slab_geometry()
    rng_b = np.random.default_rng(seed)
    pos_b, dirv_b, energy_b = _make_beam(n, energy_keV)
    rec = TrajectoryRecorder()
    result_recorded = transport_photons(pos_b, dirv_b, energy_b, geom_b, rng_b, recorder=rec)

    assert np.array_equal(result_plain.absorbed, result_recorded.absorbed)
    assert np.array_equal(result_plain.escaped, result_recorded.escaped)
    assert np.array_equal(result_plain.n_scatter, result_recorded.n_scatter)
    assert np.array_equal(result_plain.final_energy, result_recorded.final_energy)
    assert np.allclose(pos_a, pos_b)
    assert np.allclose(dirv_a, dirv_b)
    assert rec.photon_ids


def test_trace_html_contains_trajectories_and_checkbox(tmp_path):
    """HTMLスモーク: trace出力にtrajectories/ckTrajが入り、previewは空配列になる。"""
    scene = validate_scene({"source": _BASE_SOURCE, "geometry": _BASE_GEOMETRY})
    assert scene.ok

    geom = Geometry(scene.raw["geometry"])
    rng = np.random.default_rng(5)
    pos, dirv, energy = sample_source_photons(scene.raw["source"], 20, rng)
    rec = TrajectoryRecorder()
    transport_photons(pos, dirv, energy, geom, rng, recorder=rec)
    trajectories = trajectories_to_json(rec)
    assert trajectories

    trace_out = tmp_path / "trace.html"
    write_html(scene, str(trace_out), trajectories=trajectories)
    trace_html = trace_out.read_text(encoding="utf-8")
    assert '"trajectories"' in trace_html
    assert 'ckTraj' in trace_html
    assert '"points"' in trace_html

    preview_out = tmp_path / "preview.html"
    write_html(scene, str(preview_out))
    preview_html = preview_out.read_text(encoding="utf-8")
    assert '"trajectories": []' in preview_html
    assert 'ckTraj' in preview_html  # チェックボックス自体は常に存在（空なら非表示にするだけ）
