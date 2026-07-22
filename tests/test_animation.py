"""光子軌跡アニメーション向けデータ整備（Phase 0）のテスト。

docs/plan_photon_animation.md Phase 0 に対応。
classify_trajectory の全カテゴリ検証と、materials追加が既存trace出力に
影響しないことの回帰確認を行う。
"""
from __future__ import annotations

import numpy as np

from chatcarlo.geometry import Geometry
from chatcarlo.scene import validate_scene
from chatcarlo.source import sample_source_photons
from chatcarlo.trajectory import (TrajectoryRecorder, classify_trajectory,
                                   trajectories_to_json)
from chatcarlo.transport import transport_photons

_BASE_SOURCE = {
    "kvp": 80, "position": [0, -50, 0], "direction": [0, 1, 0],
    "field": {"size_cm": [10, 10], "sid_cm": 50}, "filtration_mm_al": 2.5,
}
_BASE_GEOMETRY = [{
    "name": "slab", "shape": "box", "material": "water",
    "center": [0, 0, 0], "size_cm": [10, 10, 10],
}]


def _traj(points, events):
    """手組みのtrajectories_to_json形式トラジェクトリ（materials/summaryはダミー）。"""
    return {"points": points, "events": events,
            "energies": [80.0] * len(events), "materials": ["water"] * len(events)}


def test_classify_direct_transmission():
    """無散乱でそのまま脱出: category=direct, n_interactions=0, backscattered=False。"""
    traj = _traj(
        points=[[0, 0, 0], [10, 0, 0]],
        events=["escape"],
    )
    summary = classify_trajectory(traj)
    assert summary == {"n_interactions": 0, "fate": "escaped",
                        "backscattered": False, "category": "direct"}


def test_classify_scattered_escape():
    """散乱(コンプトン)後、前方寄りの方向のまま脱出: category=scattered_escape。"""
    traj = _traj(
        points=[[0, 0, 0], [5, 0, 0], [5, 5, 0]],
        events=["compton", "escape"],
    )
    summary = classify_trajectory(traj)
    assert summary["n_interactions"] == 1
    assert summary["fate"] == "escaped"
    assert summary["backscattered"] is False
    assert summary["category"] == "scattered_escape"


def test_classify_absorbed_after_multiple_scatters():
    """コンプトン2回の後、光電吸収: category=absorbed。"""
    traj = _traj(
        points=[[0, 0, 0], [5, 0, 0], [5, 3, 0], [5, 3, 1]],
        events=["compton", "compton", "photoelectric"],
    )
    summary = classify_trajectory(traj)
    assert summary["n_interactions"] == 3
    assert summary["fate"] == "absorbed"
    assert summary["category"] == "absorbed"


def test_classify_backscatter():
    """散乱後に入射方向へ戻って脱出: category=backscatter, backscattered=True。"""
    traj = _traj(
        points=[[0, 0, 0], [5, 0, 0], [-5, 0, 0]],
        events=["compton", "escape"],
    )
    summary = classify_trajectory(traj)
    assert summary["backscattered"] is True
    assert summary["category"] == "backscatter"


def test_classify_fluorescence_counts_as_interaction():
    """蛍光再放出(fluorescence)はn_interactionsに数える（散乱ではないが反応区間）。"""
    traj = _traj(
        points=[[0, 0, 0], [5, 0, 0], [5, 0, 1], [5, 0, 2]],
        events=["fluorescence", "compton", "escape"],
    )
    summary = classify_trajectory(traj)
    assert summary["n_interactions"] == 2
    assert summary["fate"] == "escaped"


def _run_trace(n, seed):
    scene = validate_scene({"source": _BASE_SOURCE, "geometry": _BASE_GEOMETRY})
    assert scene.ok
    geom = Geometry(scene.raw["geometry"])
    rng = np.random.default_rng(seed)
    pos, dirv, energy = sample_source_photons(scene.raw["source"], n, rng)
    rec = TrajectoryRecorder()
    transport_photons(pos, dirv, energy, geom, rng, recorder=rec)
    return trajectories_to_json(rec)


def test_materials_key_present_and_aligned_with_energies():
    """materialsはenergiesと同長で、各要素が空文字でない材料名を持つ。"""
    trajs = _run_trace(n=50, seed=42)
    assert trajs
    for t in trajs:
        assert len(t["materials"]) == len(t["energies"])
        assert all(m for m in t["materials"])


def test_summary_key_present_for_every_trajectory():
    """全トラジェクトリにsummaryが付与され、既知のcategory値のみを取る。"""
    trajs = _run_trace(n=50, seed=42)
    assert trajs
    valid_categories = {"direct", "scattered_escape", "backscatter", "absorbed"}
    for t in trajs:
        assert "summary" in t
        assert t["summary"]["category"] in valid_categories
        assert t["summary"]["fate"] in ("absorbed", "escaped")


def test_animate_write_html_smoke(tmp_path):
    """chatcarlo.animate: HTML生成のスモークテスト（Phase 1）。

    転帰カテゴリ集計・再生用データ埋め込み・外部URL不参照を確認する。
    """
    from chatcarlo.animate import categorize_trajectories, write_html

    trajs = _run_trace(n=100, seed=42)
    counts = categorize_trajectories(trajs)
    assert sum(counts.values()) == len(trajs)
    assert set(counts) == {"direct", "scattered_escape", "backscatter", "absorbed"}

    scene = validate_scene({"source": _BASE_SOURCE, "geometry": _BASE_GEOMETRY})
    out = tmp_path / "anim.html"
    write_html(scene, str(out), trajs, title="test animation")
    html = out.read_text(encoding="utf-8")
    assert '"trajectories"' in html
    assert '"summary"' in html
    assert "buildFrames" in html
    assert "buildLegs" in html and "projectCam" in html
    assert "qSlerp" in html and "computeFirstPersonCam" in html
    assert "http://" not in html and "https://" not in html


def test_animate_html_phase4_ui_present(tmp_path):
    """chatcarlo.animate: Phase 4（演出仕上げ）のUI要素がHTMLに存在する。"""
    from chatcarlo.animate import write_html

    trajs = _run_trace(n=50, seed=42)
    scene = validate_scene({"source": _BASE_SOURCE, "geometry": _BASE_GEOMETRY})
    out = tmp_path / "anim_phase4.html"
    write_html(scene, str(out), trajs, title="test animation")
    html = out.read_text(encoding="utf-8")

    assert 'id="chkRaw"' in html and 'id="rawPanel"' in html
    assert 'id="chkSound"' in html and "playBlip" in html and "AudioContext" in html
    assert 'id="chkColor"' in html and "colorEnabled" in html
    assert 'id="rngSlowmo"' in html and 'id="rngPause"' in html
    assert "理解のための演出です" in html


def test_raw_data_panel_values_match_recorder_log(tmp_path):
    """生データテーブルは演出を経由せず、recorderログの値をそのまま埋め込む
    （区間ごとのmaterials/energies/eventsは元データと完全一致、区間長・散乱角は
    座標から導出されるが物理量自体は再計算しない）。
    """
    import json
    import re

    from chatcarlo.animate import write_html

    trajs = _run_trace(n=50, seed=42)
    scene = validate_scene({"source": _BASE_SOURCE, "geometry": _BASE_GEOMETRY})
    out = tmp_path / "anim_rawcheck.html"
    write_html(scene, str(out), trajs, title="test animation")
    html = out.read_text(encoding="utf-8")

    m = re.search(r"const DATA = (.*);\n", html)
    embedded = json.loads(m.group(1))["trajectories"]
    assert len(embedded) == len(trajs)
    for orig, emb in zip(trajs, embedded):
        assert emb["materials"] == orig["materials"]
        assert emb["energies"] == orig["energies"]
        assert emb["events"] == orig["events"]
        assert emb["points"] == orig["points"]


def test_recorder_extension_does_not_change_core_polyline_regression():
    """回帰確認: material追加後もpoints/energies/eventsの中身は既存仕様どおり
    （区間の連続性・終端イベント制約）を満たし続ける。materials/summaryは
    付加情報であり、既存キーの値には影響しない。
    """
    trajs = _run_trace(n=200, seed=42)
    assert trajs
    for t in trajs:
        events = t["events"]
        assert events[-1] in ("photoelectric", "escape")
        assert "photoelectric" not in events[:-1]
        assert len(t["points"]) == len(t["energies"]) + 1


def test_shared_canvas_js_is_a_single_source_not_duplicated(tmp_path):
    """preview.pyとanimate.pyの俯瞰カメラ描画コード（projectOrbit/poly/seg/
    drawAxes）はPython側の単一の文字列(_SHARED_CANVAS_JS)から両テンプレート
    へ注入される。手作業でのJSコピペ複製が再発していないか、生成後の
    HTMLに各関数がちょうど1回だけ現れることを確認する（かつてstr.replace()が
    コメント内の"__SHARED_CANVAS_JS__"という文言まで置換してしまい、
    二重に埋め込まれるバグがあった）。
    """
    import chatcarlo.animate as animate_mod
    import chatcarlo.preview as preview_mod
    from chatcarlo.animate import write_html as write_animate_html
    from chatcarlo.preview import write_html as write_preview_html

    assert animate_mod._SHARED_CANVAS_JS is preview_mod._SHARED_CANVAS_JS

    scene = validate_scene({"source": _BASE_SOURCE, "geometry": _BASE_GEOMETRY})
    shared_fn_markers = ("function projectOrbit(", "function poly(",
                         "function seg(", "function drawAxes(")

    preview_out = tmp_path / "preview.html"
    write_preview_html(scene, str(preview_out))
    preview_html = preview_out.read_text(encoding="utf-8")
    for marker in shared_fn_markers:
        assert preview_html.count(marker) == 1, marker
    assert "__SHARED_CANVAS_JS__" not in preview_html

    trajs = _run_trace(n=20, seed=42)
    anim_out = tmp_path / "anim.html"
    write_animate_html(scene, str(anim_out), trajs)
    anim_html = anim_out.read_text(encoding="utf-8")
    for marker in shared_fn_markers:
        assert anim_html.count(marker) == 1, marker
    assert "__SHARED_CANVAS_JS__" not in anim_html
