"""chatcarlo.animate が生成するJS再生ロジックの挙動テスト。

tests/test_animation.py のPhase 1-4テストは「該当する関数名が文字列として
HTMLに含まれるか」しか見ておらず、buildFrames/buildLegs/四元数slerp/
一人称カメラの実際の計算結果は検証していなかった（セッションの振り返りで
指摘された穴）。ここでは実際に出荷されるJSソースから対象関数をそのまま
抽出し、node で実行して数値的に正しいことを確認する。

抽出したソースを直接実行するため、関数の実装を変更すればこのテストも
追従する（コピー実装をテスト側で保持する二重管理を避ける）。
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess

import pytest

from chatcarlo.animate import render_html
from chatcarlo.scene import validate_scene

pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="node が見つからない")

_BASE_SOURCE = {
    "kvp": 80, "position": [0, -50, 0], "direction": [0, 1, 0],
    "field": {"size_cm": [10, 10], "sid_cm": 50}, "filtration_mm_al": 2.5,
}
_BASE_GEOMETRY = [{
    "name": "slab", "shape": "box", "material": "water",
    "center": [0, 0, 0], "size_cm": [10, 10, 10],
}]

# 依存順に注意する必要はない（_extract_symbolsが出現位置でソートする）。
_SYMBOLS = [
    "vsub", "vadd", "vdot", "vcross", "vnorm", "vnormalize", "vscale",
    "FP_REF_AXIS", "qFromVectors", "qFromDir", "qSlerp", "qRotate",
    "transportUp", "IDENTITY_Q", "buildUpAtSegStart", "distance", "lerp3", "TUNE",
    "INTERACTION_EVENTS", "buildFrames", "turnForward", "buildLegs",
    "legForSeg", "buildSegStarts", "computeFirstPersonCam",
    "RELAY_FOV", "RELAY_NEAR", "buildBasis", "projectCam",
]


def _find_function(src: str, name: str):
    m = re.search(rf"function {re.escape(name)}\(", src)
    if not m:
        return None
    start = m.start()
    brace = src.index("{", start)
    depth, i = 1, brace + 1
    while depth > 0:
        if src[i] == "{":
            depth += 1
        elif src[i] == "}":
            depth -= 1
        i += 1
    return start, i


def _find_const(src: str, name: str):
    m = re.search(rf"\bconst {re.escape(name)}\b", src)
    if not m:
        return None
    start = m.start()
    depth, i = 0, start
    while True:
        c = src[i]
        if c in "{[(":
            depth += 1
        elif c in "}])":
            depth -= 1
        elif c == ";" and depth == 0:
            i += 1
            break
        i += 1
    return start, i


def _extract_symbols(src: str, names: list[str]) -> str:
    """出荷JSから指定シンボル(関数/const宣言)をそのまま抜き出し、
    ソース中の出現順に並べ直して結合する（前方参照エラーを避けるため）。
    """
    spans = []
    for name in names:
        span = _find_function(src, name) or _find_const(src, name)
        assert span is not None, f"symbol not found in shipped JS: {name}"
        spans.append(span)
    spans.sort(key=lambda s: s[0])
    return "\n".join(src[a:b] for a, b in spans) + "\n"


def _render_js_module() -> str:
    scene = validate_scene({"source": _BASE_SOURCE, "geometry": _BASE_GEOMETRY})
    html = render_html(scene, trajectories=[], title="js logic test")
    m = re.search(r"<script>\n(.*)</script>", html, re.S)
    return _extract_symbols(m.group(1), _SYMBOLS)


def _run_node(js_body: str) -> None:
    script = js_body + "\n" + _ASSERTIONS
    proc = subprocess.run(["node", "-e", script], capture_output=True, text=True, timeout=30)
    if proc.returncode != 0:
        pytest.fail(f"JS logic assertions failed:\n{proc.stdout}\n{proc.stderr}")


_ASSERTIONS = r"""
const cv = {clientWidth: 800, clientHeight: 450};   // buildBasis/projectCamはDOM canvasの寸法しか見ない

const failures = [];
function check(name, cond, detail) {
  if (!cond) failures.push(name + (detail !== undefined ? ': ' + JSON.stringify(detail) : ''));
}
function vecClose(a, b, tol) {
  return Math.hypot(a[0]-b[0], a[1]-b[1], a[2]-b[2]) < (tol || 1e-6);
}

// ---- buildFrames: 単一区間、終端escape ----
TUNE.segBaseT = 0.5; TUNE.segScaleT = 2.0; TUNE.segMinT = 0.6; TUNE.segMaxT = 2.5;
TUNE.slowmoFraction = 0.3; TUNE.slowmoRate = 0.3; TUNE.pauseDuration = 0.4;

const trajA = {points: [[0,0,0],[10,0,0]], events: ['escape'], energies: [80], materials: ['air']};
const rA = buildFrames(trajA);
check('A_frameCount', rA.frames.length === 2, rA.frames.map(f => f.seg));
check('A_duration', Math.abs(rA.duration - 2.9) < 1e-9, rA.duration);
check('A_frame0_bounds', rA.frames[0].t0 === 0 && Math.abs(rA.frames[0].t1 - 2.5) < 1e-9);
check('A_frame1_pause', rA.frames[1].pause === true && rA.frames[1].event === 'escape');
check('A_frame1_bounds', Math.abs(rA.frames[1].t0 - 2.5) < 1e-9 && Math.abs(rA.frames[1].t1 - 2.9) < 1e-9);
check('A_frame1_frozen', vecClose(rA.frames[1].from, [10,0,0]) && vecClose(rA.frames[1].to, [10,0,0]));

// ---- buildFrames: 相互作用区間のスローモー分割比率 ----
const trajB = {points: [[0,0,0],[10,0,0],[10,10,0]], events: ['compton','escape'],
               energies: [80,60], materials: ['air','air']};
const rB = buildFrames(trajB);
const seg0 = rB.frames.filter(f => f.seg === 0);
check('B_seg0_count', seg0.length === 3, seg0.length);
check('B_seg0_normalDur', Math.abs((seg0[0].t1-seg0[0].t0) - 1.75) < 1e-9, seg0[0].t1-seg0[0].t0);
check('B_seg0_slowDur', Math.abs((seg0[1].t1-seg0[1].t0) - 2.5) < 1e-9, seg0[1].t1-seg0[1].t0);
check('B_seg0_pauseDur', Math.abs((seg0[2].t1-seg0[2].t0) - 0.4) < 1e-9);
check('B_seg0_pauseEvent', seg0[2].event === 'compton');
check('B_seg0_splitPoint', vecClose(seg0[0].to, [7,0,0]) && vecClose(seg0[1].from, [7,0,0]),
      [seg0[0].to, seg0[1].from]);
check('B_seg0_pauseFrozen', vecClose(seg0[2].from, [10,0,0]) && vecClose(seg0[2].to, [10,0,0]));
const seg1 = rB.frames.filter(f => f.seg === 1);
check('B_seg1_count', seg1.length === 2, seg1.length);

// ---- buildLegs: 境界のみ併合・反応で分割・camposの算出 ----
const trajC = {points: [[0,0,0],[1,0,0],[2,0,0],[2,1,0],[2,2,0]],
               events: ['boundary','boundary','compton','photoelectric'],
               energies: [80,80,80,50], materials: ['air','air','air','water']};
const legsC = buildLegs(trajC);
check('C_legCount', legsC.length === 3, legsC.map(l => [l.segStart, l.segEnd]));
check('C_leg0_span', legsC[0].segStart === 0 && legsC[0].segEnd === 1);
check('C_leg1_span', legsC[1].segStart === 2 && legsC[1].segEnd === 2);
check('C_leg2_span', legsC[2].segStart === 3 && legsC[2].segEnd === 3);
check('C_leg0_campos', vecClose(legsC[0].campos, [7,0,0]), legsC[0].campos);
check('C_leg1_campos', vecClose(legsC[1].campos, [7.6569,-5.6569,0], 1e-3), legsC[1].campos);
check('C_leg2_campos', vecClose(legsC[2].campos, [2,1,8], 1e-3), legsC[2].campos);

// ---- 対向・リレー視点の核心的な保証: 入射・出射の両軌跡が、直前レッグの
// 長さやこのレッグの進行度によらず、常にカメラ前方(z>0)に入り続けること。
// これはセッションの振り返りで見つかった不具合（旧設計では光子が相互作用点
// から離れるにつれ視線が振れて入射側が視野外に出ていた）に対する再発防止。
function checkRelayVisibility(angleDeg, incomingLen) {
  const rad = angleDeg * Math.PI / 180;
  const dirOut = [Math.cos(rad), Math.sin(rad), 0];
  // 直前レッグ(長さincomingLen, 方向(1,0,0))→P=[incomingLen,0,0]→このレッグ(方向dirOut, 長さ10)
  const traj = {
    points: [[0,0,0], [incomingLen,0,0], [incomingLen+10*dirOut[0], 10*dirOut[1], 0]],
    events: ['compton', 'escape'],
  };
  const legs = buildLegs(traj);
  const leg = legForSeg(legs, 1);
  const cam = buildBasis(leg.campos, leg.aimAt);
  const P = traj.points[1];
  // 入射側: Pからさらに入射方向を延長した先（直前レッグのもっと手前）
  const farIncoming = [P[0] - 1*incomingLen, P[1], P[2]];
  const farOutgoing = [P[0] + dirOut[0]*10, P[1] + dirOut[1]*10, P[2]];
  const prIn = projectCam(farIncoming, cam);
  const prOut = projectCam(farOutgoing, cam);
  check(`relay_visible_${angleDeg}deg_len${incomingLen}_in`, prIn !== null, {angleDeg, incomingLen});
  check(`relay_visible_${angleDeg}deg_len${incomingLen}_out`, prOut !== null, {angleDeg, incomingLen});
}
for (const angleDeg of [0.5, 5, 30, 90, 150, 179]) {
  for (const incomingLen of [1, 10, 20, 200]) {
    checkRelayVisibility(angleDeg, incomingLen);
  }
}

// ---- 四元数: 往復・slerp端点・対蹠点での非NaN ----
const dirs = [[1,0,0],[0,1,0],[0,0,1],[0,0,-1],[0.6,0.8,0],[-1,0,0],[0.267,0.535,0.802]];
for (const d of dirs) {
  const dn = vnormalize(d);
  const back = vnormalize(qRotate(qFromDir(dn), FP_REF_AXIS));
  check('Q_roundtrip_' + dn.join(','), vecClose(back, dn, 1e-6), back);
}
const qa = qFromDir([1,0,0]), qb = qFromDir([0,1,0]);
check('Q_slerp_t0', vecClose(qRotate(qSlerp(qa,qb,0), FP_REF_AXIS), [1,0,0]));
check('Q_slerp_t1', vecClose(qRotate(qSlerp(qa,qb,1), FP_REF_AXIS), [0,1,0]));
const qMid = vnormalize(qRotate(qSlerp(qa,qb,0.5), FP_REF_AXIS));
check('Q_slerp_mid_symmetry', Math.abs(vdot(qMid,[1,0,0]) - vdot(qMid,[0,1,0])) < 1e-6, qMid);
check('Q_slerp_mid_unit', Math.abs(vnorm(qMid) - 1) < 1e-6);

// ---- transportUp: 直交性と縮退フォールバック ----
const up1 = transportUp([0,0,1], [1,0,0]);
check('U_orthogonal', Math.abs(vdot(up1, [1,0,0])) < 1e-9, vdot(up1,[1,0,0]));
check('U_unit', Math.abs(vnorm(up1) - 1) < 1e-9);
const upDeg = transportUp([0,0,1], [0,0,1]);
check('U_degenerate_finite', isFinite(upDeg[0]) && isFinite(upDeg[1]) && isFinite(upDeg[2]), upDeg);
check('U_degenerate_orth', Math.abs(vdot(upDeg, [0,0,1])) < 1e-6, upDeg);
check('U_degenerate_unit', Math.abs(vnorm(upDeg) - 1) < 1e-6);

// ---- computeFirstPersonCam: 転換前後の向き・slerp境界・vignette ----
TUNE.slerpDuration = 0.3; TUNE.slerpDurationBig = 0.5; TUNE.fpBackOffset = 2.0;
const trajD = {points: [[0,0,0],[10,0,0],[10,10,0],[20,10,0]],
               events: ['compton','escape'], energies: [80,60], materials: ['air','air']};
const framesD = buildFrames(trajD).frames;
const segStartsD = buildSegStarts(framesD);
const upAtSegStartD = buildUpAtSegStart(trajD);

const d0 = computeFirstPersonCam(trajD, segStartsD, upAtSegStartD, 0.0, 0, [0,0,0]);
check('FP_seg0_forward', vecClose(d0.cam.forward, [1,0,0]), d0.cam.forward);
check('FP_seg0_campos', vecClose(d0.cam.pos, [-2,0,0]), d0.cam.pos);
check('FP_seg0_up_orthogonal', Math.abs(vdot(d0.cam.up, d0.cam.forward)) < 1e-9, d0.cam.up);

const t1start = segStartsD[1];
const d1a = computeFirstPersonCam(trajD, segStartsD, upAtSegStartD, t1start, 1, [10,0,0]);
check('FP_seg1_start_forward_eq_prevDir', vecClose(d1a.cam.forward, [1,0,0], 1e-6), d1a.cam.forward);
check('FP_seg1_start_vignette_zero', d1a.vignette === 0, d1a.vignette);

// 直角(90度)転換なのでangleDeg>90はfalse -> slerpDuration(0.3s)を使う
const d1b = computeFirstPersonCam(trajD, segStartsD, upAtSegStartD, t1start + 0.31, 1, [10,5,0]);
check('FP_seg1_end_forward_eq_nextDir', vecClose(d1b.cam.forward, [0,1,0], 1e-3), d1b.cam.forward);
check('FP_seg1_end_up_orthogonal', Math.abs(vdot(d1b.cam.up, d1b.cam.forward)) < 1e-9, d1b.cam.up);

// 大角度(>90度)転換: slerpDurationBig(0.5s)とvignette>0を確認
const trajE = {points: [[0,0,0],[10,0,0],[4,8,0]], events: ['compton','escape'],
               energies: [80,60], materials: ['air','air']};
const framesE = buildFrames(trajE).frames;
const segStartsE = buildSegStarts(framesE);
const upAtSegStartE = buildUpAtSegStart(trajE);
const eMid = computeFirstPersonCam(trajE, segStartsE, upAtSegStartE, segStartsE[1] + 0.25, 1, [7,4,0]);
check('FP_bigangle_vignette_positive', eMid.vignette > 0.3, eMid.vignette);
const eEnd = computeFirstPersonCam(trajE, segStartsE, upAtSegStartE, segStartsE[1] + 0.51, 1, [4,8,0]);
check('FP_bigangle_end_forward', vecClose(eEnd.cam.forward, vnormalize([-6,8,0]), 1e-3), eEnd.cam.forward);

// ---- 冪等性: draw()の呼び出し順序やシークバーのスクラブ方向に関わらず、
// 同じclockに対しては常に同じ姿勢が返ること（旧fpStateはdraw()履歴に
// 依存する可変状態だったため、逆再生すると順再生と異なるロールになる
// 不具合があった——その再発防止テスト）。3区間の方向を互いに直交させて
// あるのが重要: 同じ方向への射影は複数回適用しても結果が変わらず
// （冪等）、蓄積型の実装でも偶然一致してテストをすり抜けてしまうため
// （実際に単純な蓄積バグで再現し、これに気づいて区間方向を直交にした）。
const trajF = {points: [[0,0,0],[10,0,0],[10,10,0],[10,10,10]],
               events: ['compton', 'compton', 'escape'],
               energies: [80,60,50], materials: ['air','air','air']};
const framesF = buildFrames(trajF).frames;
const segStartsF = buildSegStarts(framesF);
const upAtSegStartF = buildUpAtSegStart(trajF);
const t1F = segStartsF[1], t2F = segStartsF[2];

function segAndPosF(c) {
  if (c < t1F) return [0, [c*10/t1F, 0, 0]];
  if (c < t2F) return [1, [10, (c-t1F)*10/(t2F-t1F), 0]];
  return [2, [10, 10, 0]];
}
const probeClocks = [0.05, t1F + 0.05, t1F + 0.15, t2F + 0.05, t2F + 0.31];
const forwardOrder = probeClocks.map(c => {
  const [seg, pos] = segAndPosF(c);
  return computeFirstPersonCam(trajF, segStartsF, upAtSegStartF, c, seg, pos);
});
const shuffled = [2, 0, 4, 1, 3];   // 時系列順ではない呼び出し順（逆再生・ジャンプを模す）
for (const idx of shuffled) {
  const c = probeClocks[idx];
  const [seg, pos] = segAndPosF(c);
  const r = computeFirstPersonCam(trajF, segStartsF, upAtSegStartF, c, seg, pos);
  check(`idempotent_clock${idx}`,
        vecClose(r.cam.up, forwardOrder[idx].cam.up, 1e-12) &&
        vecClose(r.cam.forward, forwardOrder[idx].cam.forward, 1e-12),
        {shuffledUp: r.cam.up, expectedUp: forwardOrder[idx].cam.up});
}

if (failures.length) {
  console.log(JSON.stringify(failures, null, 2));
  process.exit(1);
} else {
  console.log('ALL_PASS');
}
"""


def test_js_playback_and_camera_logic_matches_expected_math():
    """出荷JSのbuildFrames/buildLegs/四元数/一人称カメラを、実データで数値検証する。"""
    module_src = _render_js_module()
    _run_node(module_src)
