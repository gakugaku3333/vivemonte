"""教育用の光子軌跡アニメーション — scene.yaml + 軌跡ログから自己完結HTMLを生成する。

docs/plan_photon_animation.md Phase 1〜4に対応（俯瞰/対向リレー/一人称の3カメラ、
演出強度設定、生データ表示モード）。`chatcarlo/preview.py` と同じくCDN非依存の
vanilla JS + canvasで実装する。ジオメトリー・ビームのメッシュ化とJSON整形は
preview.scene_to_json をそのまま再利用し、このモジュールが足すのは
「時間軸に沿った再生」だけ。

物理量（位置・エネルギー・相互作用種別）はJS側で一切再計算しない。演出
（時間の圧縮・スローモー・色）とデータの境界を守るのは docs/plan_photon_animation.md
§1.4 の禁止事項。
"""
from __future__ import annotations

import json

from .preview import _SHARED_CANVAS_JS, scene_to_json
from .scene import Scene

CATEGORY_LABELS = {
    "direct": "直接透過",
    "scattered_escape": "散乱後脱出",
    "backscatter": "後方散乱",
    "absorbed": "吸収",
}
CATEGORY_ORDER = ["direct", "scattered_escape", "backscatter", "absorbed"]


def categorize_trajectories(trajectories: list[dict]) -> dict[str, int]:
    """転帰カテゴリ別の件数集計（CLI表示・対話AI導線の情報源）。"""
    counts = {k: 0 for k in CATEGORY_ORDER}
    for t in trajectories:
        cat = t["summary"]["category"]
        counts[cat] = counts.get(cat, 0) + 1
    return counts


def render_html(scene: Scene, trajectories: list[dict],
                 title: str = "ChatCarlo photon animation") -> str:
    data = scene_to_json(scene, trajectories=trajectories)
    # カテゴリ名→ラベルの対応はcategorize_trajectories/CLI表示と共通の
    # CATEGORY_LABELS/CATEGORY_ORDERから注入する（3箇所に手で複製しない）。
    return (_TEMPLATE.replace("__TITLE__", title)
            .replace("__DATA__", json.dumps(data, ensure_ascii=False))
            .replace("__CATEGORY_LABELS__", json.dumps(CATEGORY_LABELS, ensure_ascii=False))
            .replace("__CATEGORY_ORDER__", json.dumps(CATEGORY_ORDER, ensure_ascii=False))
            .replace("__SHARED_CANVAS_JS__", _SHARED_CANVAS_JS))


def write_html(scene: Scene, out_path: str, trajectories: list[dict],
               title: str = "ChatCarlo photon animation") -> str:
    html = render_html(scene, trajectories, title=title)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    return out_path


_TEMPLATE = r"""<!doctype html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  :root { --bg:#0f172a; --panel:#1e293b; --text:#e2e8f0; --dim:#94a3b8; --accent:#fbbf24; }
  @media (prefers-color-scheme: light) {
    :root { --bg:#f1f5f9; --panel:#ffffff; --text:#0f172a; --dim:#64748b; }
  }
  * { margin:0; box-sizing:border-box; }
  body { background:var(--bg); color:var(--text); font-family:system-ui,-apple-system,sans-serif;
         display:flex; flex-direction:column; height:100vh; overflow:hidden; }
  header { padding:8px 14px; display:flex; gap:10px; align-items:center; flex-wrap:wrap;
           border-bottom:1px solid #47556933; }
  header h1 { font-size:14px; font-weight:600; margin-right:auto; }
  button { background:var(--panel); color:var(--text); border:1px solid #47556933;
           border-radius:6px; padding:5px 11px; font-size:12px; cursor:pointer; }
  button:hover { border-color:var(--accent); }
  button.active { border-color:var(--accent); color:var(--accent); }
  button:disabled { opacity:.4; cursor:default; }
  #main { flex:1; display:flex; min-height:0; }
  #side { width:270px; overflow-y:auto; border-right:1px solid #47556933; padding:10px; font-size:12px; }
  #side h2 { font-size:11px; color:var(--dim); margin:10px 0 4px; text-transform:uppercase; letter-spacing:.04em; }
  #side h2:first-child { margin-top:0; }
  .photon-item { padding:5px 8px; border-radius:5px; cursor:pointer; display:flex;
                 justify-content:space-between; gap:6px; }
  .photon-item:hover { background:#47556933; }
  .photon-item.selected { background:var(--accent); color:#1a1300; }
  .photon-item .meta { color:var(--dim); font-size:10.5px; }
  .photon-item.selected .meta { color:#1a1300cc; }
  #wrap { flex:1; position:relative; }
  canvas { position:absolute; inset:0; width:100%; height:100%; cursor:grab; }
  #hud { position:absolute; top:10px; left:12px; background:var(--panel); border-radius:8px;
         padding:10px 14px; font-size:12px; line-height:1.9; opacity:.95; min-width:190px; }
  #hud .row { display:flex; justify-content:space-between; gap:16px; }
  #hud .label { color:var(--dim); }
  #playbar { position:absolute; bottom:0; left:0; right:0; background:var(--panel);
             padding:8px 14px; display:flex; align-items:center; gap:10px; font-size:12px; }
  #playbar input[type=range] { flex:1; }
  #hint { position:absolute; bottom:44px; left:12px; font-size:11px; color:var(--dim); }
  #note { position:absolute; top:10px; right:12px; max-width:250px; background:#1e293bcc;
          border:1px solid #47556955; border-radius:8px; padding:8px 10px; font-size:10.5px;
          line-height:1.5; color:var(--dim); }
  #placeholder { position:absolute; inset:0; display:flex; align-items:center; justify-content:center;
                 color:var(--dim); font-size:13px; }
  #settingsPanel { position:absolute; top:44px; right:12px; background:var(--panel);
                   border:1px solid #47556955; border-radius:8px; padding:10px 12px;
                   font-size:11.5px; display:none; min-width:210px; z-index:2; }
  #settingsPanel label { display:flex; align-items:center; gap:6px; cursor:pointer;
                         margin:6px 0; color:var(--text); }
  #settingsPanel .sliderRow { display:flex; flex-direction:column; gap:2px; margin:8px 0; }
  #settingsPanel .sliderRow span { color:var(--dim); font-size:10.5px; }
  #settingsPanel input[type=range] { width:100%; }
  #rawPanel { position:absolute; bottom:52px; right:12px; width:300px; max-height:60%;
              overflow-y:auto; background:var(--panel); border-radius:8px; padding:8px 10px;
              font-size:10.5px; display:none; }
  #rawPanel table { border-collapse:collapse; width:100%; }
  #rawPanel th, #rawPanel td { padding:2px 5px; text-align:right; white-space:nowrap; }
  #rawPanel th:first-child, #rawPanel td:first-child { text-align:left; }
  #rawPanel th { color:var(--dim); font-weight:600; }
  #rawPanel tr.current { background:var(--accent); color:#1a1300; }
</style>
</head>
<body>
<header>
  <h1>__TITLE__</h1>
  <span id="camModeGroup">
    <button data-cam="orbit" class="active">俯瞰</button>
    <button data-cam="relay">対向リレー</button>
    <button data-cam="fp">一人称</button>
  </span>
  <span id="orbitViewGroup">
    <button data-view="iso">等角</button>
    <button data-view="front">正面 (−Y)</button>
    <button data-view="side">側面 (+X)</button>
    <button data-view="top">上面 (+Z)</button>
  </span>
  <button id="btnSettings">⚙ 設定</button>
</header>
<div id="main">
  <div id="side"></div>
  <div id="wrap">
    <canvas id="cv"></canvas>
    <div id="hud" style="display:none"></div>
    <div id="note">⚠ 時間・色・カメラの動きは理解のための演出です<br>（軌跡・エネルギー・相互作用は計算ログの実データ）</div>
    <div id="hint">ドラッグ: 回転　/　Shift+ドラッグ: 平行移動　/　ホイール: ズーム</div>
    <div id="placeholder">左のリストから光子を選んでください</div>
    <div id="settingsPanel">
      <label><input type="checkbox" id="chkRaw"> 生データテーブルを表示</label>
      <label><input type="checkbox" id="chkSound"> 効果音（相互作用時にピッチ変化）</label>
      <label><input type="checkbox" id="chkColor" checked> エネルギーで色分け</label>
      <div class="sliderRow">
        <span>スローモー倍率: <span id="lblSlowmo">0.30</span>×</span>
        <input type="range" id="rngSlowmo" min="0.1" max="1.0" step="0.05" value="0.3">
      </div>
      <div class="sliderRow">
        <span>相互作用の停止時間: <span id="lblPause">0.40</span>s</span>
        <input type="range" id="rngPause" min="0.1" max="1.0" step="0.05" value="0.4">
      </div>
    </div>
    <div id="rawPanel"></div>
    <div id="playbar" style="display:none">
      <button id="btnPlay">▶</button>
      <button id="btnRestart">⟲</button>
      <input type="range" id="seek" min="0" max="1000" value="0">
      <span id="speedGroup">
        <button data-speed="0.5">0.5×</button>
        <button data-speed="1" class="active">1×</button>
        <button data-speed="2">2×</button>
      </span>
    </div>
  </div>
</div>
<script>
const DATA = __DATA__;
const cv = document.getElementById('cv'), ctx = cv.getContext('2d');
let yaw = -0.7, pitch = 0.42, dist = DATA.radius * 2.6;
let panX = 0, panY = 0, dragging = false, panning = false, lx = 0, ly = 0;

const views = { iso:[-0.7,0.42], front:[0, 0], side:[Math.PI/2, 0], top:[0, Math.PI/2 - 1e-3] };
document.querySelectorAll('button[data-view]').forEach(b => b.onclick = () => {
  [yaw, pitch] = views[b.dataset.view]; panX = panY = 0; draw();
});

// ---- カメラモード（俯瞰 / 対向リレー） ----
let cameraMode = 'orbit';
document.querySelectorAll('#camModeGroup button').forEach(b => b.onclick = () => {
  document.querySelectorAll('#camModeGroup button').forEach(x => x.classList.remove('active'));
  b.classList.add('active');
  cameraMode = b.dataset.cam;
  document.getElementById('orbitViewGroup').style.display = cameraMode === 'orbit' ? '' : 'none';
  document.getElementById('hint').style.display = cameraMode === 'orbit' ? '' : 'none';
  draw();
});

// ---- 演出強度設定パネル ----
let colorEnabled = true, soundEnabled = false, audioCtx = null;
document.getElementById('btnSettings').onclick = () => {
  const p = document.getElementById('settingsPanel');
  p.style.display = p.style.display === 'block' ? 'none' : 'block';
};
document.getElementById('chkColor').onchange = e => { colorEnabled = e.target.checked; draw(); };
document.getElementById('chkSound').onchange = e => {
  soundEnabled = e.target.checked;
  if (soundEnabled && !audioCtx) {
    audioCtx = new (window.AudioContext || window.webkitAudioContext)();
  }
};
document.getElementById('chkRaw').onchange = e => {
  document.getElementById('rawPanel').style.display = e.target.checked ? 'block' : 'none';
};
document.getElementById('rngSlowmo').oninput = e => {
  TUNE.slowmoRate = parseFloat(e.target.value);
  document.getElementById('lblSlowmo').textContent = TUNE.slowmoRate.toFixed(2);
  rebuildCurrentTimeline();
};
document.getElementById('rngPause').oninput = e => {
  TUNE.pauseDuration = parseFloat(e.target.value);
  document.getElementById('lblPause').textContent = TUNE.pauseDuration.toFixed(2);
  rebuildCurrentTimeline();
};

// 効果音: 相互作用の一時停止フレームに入った瞬間だけ鳴らす（再生ボタン初回操作で
// AudioContextが生成されている場合のみ）。ピッチはエネルギーに応じて変化させる。
let lastSoundFrame = null;
function playBlip(energyKeV) {
  if (!soundEnabled || !audioCtx) return;
  const freq = 200 + 1800 * (energyKeV / eMaxGlobal);
  const osc = audioCtx.createOscillator(), gain = audioCtx.createGain();
  osc.type = 'sine';
  osc.frequency.value = freq;
  gain.gain.setValueAtTime(0.15, audioCtx.currentTime);
  gain.gain.exponentialRampToValueAtTime(0.001, audioCtx.currentTime + 0.08);
  osc.connect(gain); gain.connect(audioCtx.destination);
  osc.start(); osc.stop(audioCtx.currentTime + 0.08);
}

function resize() {
  const r = cv.parentElement.getBoundingClientRect(), dpr = devicePixelRatio || 1;
  cv.width = r.width * dpr; cv.height = r.height * dpr;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0); draw();
}
addEventListener('resize', resize);

cv.onmousedown = e => { dragging = true; panning = e.shiftKey; lx = e.clientX; ly = e.clientY; cv.style.cursor='grabbing'; };
addEventListener('mouseup', () => { dragging = false; cv.style.cursor='grab'; });
addEventListener('mousemove', e => {
  if (!dragging) return;
  const dx = e.clientX - lx, dy = e.clientY - ly; lx = e.clientX; ly = e.clientY;
  if (panning) { panX += dx; panY += dy; }
  else { yaw += dx * 0.008; pitch = Math.max(-1.55, Math.min(1.55, pitch + dy * 0.008)); }
  draw();
});
cv.addEventListener('wheel', e => { e.preventDefault(); dist *= Math.exp(e.deltaY * 0.001); draw(); }, {passive:false});

// projectOrbit（俯瞰カメラの投影。preview.pyと共通のJSチャンク、Python側の
// _SHARED_CANVAS_JSから注入——ソースの二重管理を避けるため、ここにJSを直接
// 書かない）はここで定義される。poly/seg/drawAxesも同様に共通チャンク側。
__SHARED_CANVAS_JS__

// ---- 位置・姿勢つき透視カメラ（対向・リレー視点用の一般化） ----
function vsub(a, b) { return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]; }
function vadd(a, b) { return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]; }
function vdot(a, b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
function vcross(a, b) {
  return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
}
function vnorm(a) { return Math.hypot(a[0], a[1], a[2]); }
function vnormalize(a) { const n = vnorm(a) || 1; return [a[0]/n, a[1]/n, a[2]/n]; }

const RELAY_FOV = 60 * Math.PI / 180;
const RELAY_NEAR = 1e-3;

// lookAtカメラの基底構築（glm::lookAt と同じ標準的な導出）。
// forwardがワールドupとほぼ平行（真上/真下を向く）場合はフォールバック軸を使う。
function buildBasis(pos, target) {
  const forward = vnormalize(vsub(target, pos));
  let upHint = [0, 0, 1];
  let right = vcross(forward, upHint);
  if (vnorm(right) < 1e-4) { upHint = [1, 0, 0]; right = vcross(forward, upHint); }
  right = vnormalize(right);
  const up = vnormalize(vcross(right, forward));
  return {pos, right, up, forward};
}

function projectCam(p, cam) {
  const v = vsub(p, cam.pos);
  const x = vdot(v, cam.right), y = vdot(v, cam.up), z = vdot(v, cam.forward);
  if (z < RELAY_NEAR) return null;
  const w = cv.clientWidth, h = cv.clientHeight;
  const f = 0.5 * Math.min(w, h) / Math.tan(RELAY_FOV / 2);
  return [w/2 + f * x / z, h/2 - f * y / z, z];
}

// 描画コード全体（poly/seg/dot/drawAxes/drawScene）が使う投影ディスパッチャ。
// activeCamがnullなら俯瞰、そうでなければそのフレームのリレー/一人称カメラを使う。
let activeCam = null;
function project(p) { return activeCam ? projectCam(p, activeCam) : projectOrbit(p); }

// ---- 四元数（一人称視点の方向転換slerp用） ----
// [x,y,z,w] 形式。fromDir/slerp/rotateの3関数のみで足りる（docs/plan_photon_animation.md Phase 3）。
function vscale(a, s) { return [a[0]*s, a[1]*s, a[2]*s]; }
const FP_REF_AXIS = [0, 0, 1];

function qFromVectors(a, b) {
  const d = Math.max(-1, Math.min(1, vdot(a, b)));
  if (d > 1 - 1e-9) return [0, 0, 0, 1];
  if (d < -1 + 1e-9) {
    let axis = vcross([1, 0, 0], a);
    if (vnorm(axis) < 1e-6) axis = vcross([0, 1, 0], a);
    axis = vnormalize(axis);
    return [axis[0], axis[1], axis[2], 0];
  }
  const axis = vnormalize(vcross(a, b));
  const angle = Math.acos(d);
  const s = Math.sin(angle / 2);
  return [axis[0]*s, axis[1]*s, axis[2]*s, Math.cos(angle / 2)];
}
function qFromDir(dir) { return qFromVectors(FP_REF_AXIS, dir); }

function qSlerp(qa, qb, t) {
  let [x0, y0, z0, w0] = qa, [x1, y1, z1, w1] = qb;
  let d = x0*x1 + y0*y1 + z0*z1 + w0*w1;
  if (d < 0) { x1 = -x1; y1 = -y1; z1 = -z1; w1 = -w1; d = -d; }
  if (d > 0.9995) {
    const rx = x0+(x1-x0)*t, ry = y0+(y1-y0)*t, rz = z0+(z1-z0)*t, rw = w0+(w1-w0)*t;
    const n = Math.hypot(rx, ry, rz, rw) || 1;
    return [rx/n, ry/n, rz/n, rw/n];
  }
  const theta0 = Math.acos(d), theta = theta0 * t;
  const sinTheta0 = Math.sin(theta0), sinTheta = Math.sin(theta);
  const s0 = Math.cos(theta) - d * sinTheta / sinTheta0;
  const s1 = sinTheta / sinTheta0;
  return [s0*x0+s1*x1, s0*y0+s1*y1, s0*z0+s1*z1, s0*w0+s1*w1];
}

function qRotate(q, v) {
  const [x, y, z, w] = q, qv = [x, y, z];
  const uv = vcross(qv, v);
  const uuv = vcross(qv, uv);
  return [v[0]+2*(w*uv[0]+uuv[0]), v[1]+2*(w*uv[1]+uuv[1]), v[2]+2*(w*uv[2]+uuv[2])];
}

// upベクトルの平行移動的な引き継ぎ: 世界z軸には固定せず、直前のupを新forwardへ
// 直交射影する。forwardとほぼ平行（真上/真下向き）ならフォールバック軸を使う。
// 単体では「新方向へのGram-Schmidt直交化」でしかなく、回転を表さない
// （回転として引き継ぐ経路はqFromVectors側、buildUpAtSegStart参照）。
function transportUp(oldUp, newForward) {
  let up = vsub(oldUp, vscale(newForward, vdot(oldUp, newForward)));
  if (vnorm(up) < 1e-4) {
    const hint = Math.abs(newForward[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    up = vsub(hint, vscale(newForward, vdot(hint, newForward)));
  }
  return vnormalize(up);
}

const IDENTITY_Q = [0, 0, 0, 1];

// 光子選択時に一度だけ計算する「各区間の始点で決着したupベクトル」の配列。
// 反応（相互作用）による方向転換は、prevDir→nextDirの回転をupにもそのまま
// 適用する（qFromVectors経由）。これにより一人称カメラの姿勢はclockだけの
// 純粋関数になる——旧実装はdraw()が呼ばれるたびにモジュール直下の可変状態
// (旧fpState.up)を書き換えて引き継いでおり、シークバーを逆再生・ジャンプ
// すると順再生と異なるロールになる、フレームレート依存の不具合があった。
// boundary（方向不変）区間はGram-Schmidt直交化のみ行う（数値誤差の補正）。
function buildUpAtSegStart(traj) {
  const n = traj.events.length;
  const ups = [];
  let prevDir = vnormalize(vsub(traj.points[1], traj.points[0]));
  let up = transportUp([0, 0, 1], prevDir);
  ups.push(up);
  for (let i = 1; i < n; i++) {
    const dir = vnormalize(vsub(traj.points[i+1], traj.points[i]));
    up = INTERACTION_EVENTS.has(traj.events[i-1])
      ? qRotate(qFromVectors(prevDir, dir), up)
      : transportUp(up, dir);
    ups.push(up);
    prevDir = dir;
  }
  return ups;
}

// 一人称カメラの姿勢を求める。upAtSegStart（buildUpAtSegStartで事前計算済み）と
// clockだけから一意に決まる純粋関数——draw()の呼び出し履歴に依存する可変状態は
// 持たない。相互作用直後(traj.events[curSeg-1]が反応)は旧方向→新方向を
// 四元数slerpで回す（大角度転換は500msへ延長）。
function computeFirstPersonCam(traj, segStarts, upAtSegStart, clock, curSeg, curPos) {
  const nextDir = vnormalize(vsub(traj.points[curSeg+1], traj.points[curSeg]));
  let interpDir = nextDir, up = upAtSegStart[curSeg], vignette = 0;
  if (curSeg > 0 && INTERACTION_EVENTS.has(traj.events[curSeg-1])) {
    const prevDir = vnormalize(vsub(traj.points[curSeg], traj.points[curSeg-1]));
    const angleDeg = Math.acos(Math.max(-1, Math.min(1, vdot(prevDir, nextDir)))) * 180 / Math.PI;
    const dur = angleDeg > 90 ? TUNE.slerpDurationBig : TUNE.slerpDuration;
    const frac = Math.min(1, Math.max(0, (clock - segStarts[curSeg]) / dur));
    if (frac < 1) {
      const q = qSlerp(qFromDir(prevDir), qFromDir(nextDir), frac);
      interpDir = vnormalize(qRotate(q, FP_REF_AXIS));
      // upも同じ区間遷移のあいだ、prevDir→nextDirの回転をfracぶんだけ
      // 部分適用して滑らかに回す（interpDir自体とは別の回転経路だが、
      // 端点(frac=0/1)では厳密に一致し、途中は後段の再直交化で吸収する）。
      const qPartial = qSlerp(IDENTITY_Q, qFromVectors(prevDir, nextDir), frac);
      up = qRotate(qPartial, upAtSegStart[curSeg-1]);
      if (angleDeg > 90) vignette = Math.sin(Math.PI * frac) * 0.35;
    }
  }
  const right = vnormalize(vcross(interpDir, up));
  up = vnormalize(vcross(right, interpDir));
  const camPos = [curPos[0] - interpDir[0]*TUNE.fpBackOffset,
                  curPos[1] - interpDir[1]*TUNE.fpBackOffset,
                  curPos[2] - interpDir[2]*TUNE.fpBackOffset];
  return {cam: {pos: camPos, right, up, forward: interpDir}, vignette};
}

// poly/seg/drawAxesはprojectOrbitと同じ共通JSチャンク（上のprojectOrbit定義
// 箇所で注入済み）にまとめて含まれる。ここではanimate.py固有のdotのみ定義する。
function dot(a, b, r, color) {
  const p = project(a); if (!p) return;
  ctx.fillStyle = color; ctx.beginPath(); ctx.arc(p[0], p[1], r, 0, 7); ctx.fill();
}

function drawScene() {
  const all = [];
  for (const o of DATA.objects)
    for (const f of o.faces) {
      let d = 0, n = 0;
      for (const p of f) { const pr = project(p); if (!pr) { d = -1; break; } d += pr[2]; n++; }
      if (d > 0) all.push([d / n, f, o.color]);
    }
  all.sort((a, b) => b[0] - a[0]);
  for (const [, f, col] of all) poly(f, col, 0.1);
  for (const o of DATA.objects)
    for (const e of o.edges) seg(e[0], e[1], o.color, 1.0);

  const b = DATA.beam;
  const nc = b.corners.length;
  const apex = nc > 4 ? b.corners.filter((c, i) => i % 4 === 0) : b.corners;
  for (const c of apex) seg(b.source, c, '#fbbf2455', 1, [5, 4]);
  for (let i = 0; i < nc; i++) seg(b.corners[i], b.corners[(i+1)%nc], '#fbbf2477', 1.5);
}

function markerAt(p, event, col) {
  const pr = project(p); if (!pr) return;
  if (event === 'photoelectric') {
    ctx.fillStyle = col; ctx.beginPath(); ctx.arc(pr[0], pr[1], 4, 0, 7); ctx.fill();
  } else if (event === 'compton') {
    ctx.strokeStyle = col; ctx.lineWidth = 1.4;
    ctx.beginPath(); ctx.arc(pr[0], pr[1], 4, 0, 7); ctx.stroke();
  } else if (event === 'rayleigh') {
    ctx.strokeStyle = col; ctx.lineWidth = 1.4;
    ctx.beginPath();
    ctx.moveTo(pr[0], pr[1]-5); ctx.lineTo(pr[0]+5, pr[1]);
    ctx.lineTo(pr[0], pr[1]+5); ctx.lineTo(pr[0]-5, pr[1]);
    ctx.closePath(); ctx.stroke();
  } else if (event === 'escape') {
    ctx.strokeStyle = col; ctx.lineWidth = 1.6;
    ctx.beginPath();
    ctx.moveTo(pr[0]-4, pr[1]-4); ctx.lineTo(pr[0]+4, pr[1]+4);
    ctx.moveTo(pr[0]-4, pr[1]+4); ctx.lineTo(pr[0]+4, pr[1]-4);
    ctx.stroke();
  } else if (event === 'fluorescence') {
    ctx.fillStyle = col; ctx.fillRect(pr[0]-4, pr[1]-4, 8, 8);
  }
}

// ---- 転帰分類パネル ----
// ラベル・並び順はPython側(CATEGORY_LABELS/CATEGORY_ORDER、animate.pyの
// render_html)から注入する。CLI表示(__main__.py)とここで手書きの複製を
// 持たないようにするため。
const CATEGORY_LABELS = __CATEGORY_LABELS__;
const CATEGORY_ORDER = __CATEGORY_ORDER__;

function buildSidePanel() {
  const side = document.getElementById('side');
  const byCategory = {};
  DATA.trajectories.forEach((t, i) => {
    const cat = t.summary.category;
    (byCategory[cat] = byCategory[cat] || []).push(i);
  });
  let html = '';
  for (const cat of CATEGORY_ORDER) {
    const idxs = byCategory[cat] || [];
    if (!idxs.length) continue;
    html += `<h2>${CATEGORY_LABELS[cat]}（${idxs.length}）</h2>`;
    for (const i of idxs) {
      const t = DATA.trajectories[i];
      const eFinal = t.energies[t.energies.length - 1].toFixed(1);
      html += `<div class="photon-item" data-idx="${i}">` +
              `<span>光子 #${i}</span>` +
              `<span class="meta">反応${t.summary.n_interactions}回 / ${eFinal}keV</span></div>`;
    }
  }
  side.innerHTML = html;
  side.querySelectorAll('.photon-item').forEach(el => {
    el.onclick = () => selectPhoton(parseInt(el.dataset.idx, 10));
  });
}

// ---- 再生タイムライン ----
const TUNE = {
  segBaseT: 0.5, segScaleT: 2.0, segMinT: 0.6, segMaxT: 2.5,
  slowmoFraction: 0.3, slowmoRate: 0.3, pauseDuration: 0.4,
  slerpDuration: 0.3, slerpDurationBig: 0.5, fpBackOffset: 2.0,
  terminalFadeDuration: 1.0,
};
const INTERACTION_EVENTS = new Set(['photoelectric', 'compton', 'rayleigh', 'fluorescence']);

function distance(a, b) {
  return Math.hypot(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}
function lerp3(a, b, f) {
  return [a[0]+(b[0]-a[0])*f, a[1]+(b[1]-a[1])*f, a[2]+(b[2]-a[2])*f];
}

function buildFrames(traj) {
  const n = traj.events.length;
  let lMax = 0;
  for (let i = 0; i < n; i++) lMax = Math.max(lMax, distance(traj.points[i], traj.points[i+1]));
  lMax = lMax || 1;

  const frames = [];
  let t = 0;
  for (let i = 0; i < n; i++) {
    const p0 = traj.points[i], p1 = traj.points[i+1];
    const L = distance(p0, p1);
    const base = Math.min(TUNE.segMaxT, Math.max(TUNE.segMinT,
      TUNE.segBaseT + TUNE.segScaleT * Math.sqrt(L / lMax)));
    const isInteraction = INTERACTION_EVENTS.has(traj.events[i]);
    const common = {seg: i, energy: traj.energies[i], material: traj.materials[i]};
    if (isInteraction) {
      const splitFrac = 1 - TUNE.slowmoFraction;
      const pMid = lerp3(p0, p1, splitFrac);
      const normalDur = base * splitFrac;
      const slowDur = base * TUNE.slowmoFraction / TUNE.slowmoRate;
      frames.push({...common, t0: t, t1: t + normalDur, from: p0, to: pMid});
      t += normalDur;
      frames.push({...common, t0: t, t1: t + slowDur, from: pMid, to: p1});
      t += slowDur;
      frames.push({...common, t0: t, t1: t + TUNE.pauseDuration, from: p1, to: p1,
                   pause: true, event: traj.events[i]});
      t += TUNE.pauseDuration;
    } else {
      frames.push({...common, t0: t, t1: t + base, from: p0, to: p1});
      t += base;
      if (traj.events[i] === 'escape') {
        frames.push({...common, t0: t, t1: t + TUNE.pauseDuration, from: p1, to: p1,
                     pause: true, event: 'escape'});
        t += TUNE.pauseDuration;
      }
    }
  }
  return {frames, duration: t};
}

// ---- 対向・リレー視点: 「レッグ」（材料境界通過のみで方向が変わらない
// 連続区間の束）ごとにカメラを1台割り当てる。境界通過だけでカメラを
// 切り替えると無意味なカットが多発するため、方向が同じ区間は併合する。
//
// カメラ姿勢はレッグ内で完全に固定し、光子の現在位置を追い続けて
// 毎フレーム視線を振り直す実装にはしない。振り直す方式だと、光子が
// 相互作用点Pから離れるにつれて視線が出射方向へ回転していき、数フレーム
// 後には入射側の軌跡がカメラの後方（視錐台の外）へ外れて見えなくなる
// ——実測で確認した不具合（入射・出射が同一フレームに同時に映るという
// Phase 2の受け入れ基準を満たしていなかった）。
//
// 固定視線の向きは forward = normalize(dirOut − dirIn) を用いる。この
// 選び方には数学的根拠があり、入射側の軌跡（P − dirIn·s, s>0）と
// 出射側の軌跡（P + dirOut·s, s>0）のカメラ深度は
//   z_in(s)  = offset − s・dot(dirIn , forward)
//   z_out(s) = offset + s・dot(dirOut, forward)
// で、forward = normalize(dirOut−dirIn) のとき dot(dirIn,forward) ≤ 0 かつ
// dot(dirOut,forward) ≥ 0 が常に成り立つため、区間の長さ s によらず
// 両方とも z > 0（カメラ前方）を保つ。つまり直前のレッグがどれだけ
// 長くても、また今のレッグをどれだけ進んでも、入射・出射の両方向が
// 同一フレームに映り続けることが構造的に保証される（ほぼ180度の
// 後方散乱では入射・出射がほぼ同一直線になり画面上で重なるが、これは
// 物理的に正しい——ほぼ同じ経路を逆走するのだから区別できなくて当然）。
function turnForward(dirIn, dirOut) {
  const diff = vsub(dirOut, dirIn);
  if (vnorm(diff) < 1e-3) {
    // ほぼ0度（直進に近い散乱）: 有意な転換が無いので任意の垂直方向を使う
    return vnormalize(vcross(dirIn, Math.abs(dirIn[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]));
  }
  return vnormalize(diff);
}

function buildLegs(traj) {
  const pts = traj.points;
  const n = traj.events.length;
  const legs = [];
  let segStart = 0;
  for (let i = 0; i < n - 1; i++) {
    const dirI = vnormalize(vsub(pts[i+1], pts[i]));
    const dirNext = vnormalize(vsub(pts[i+2], pts[i+1]));
    const sameLeg = traj.events[i] === 'boundary' && vdot(dirI, dirNext) > 1 - 1e-6;
    if (!sameLeg) { legs.push({segStart, segEnd: i}); segStart = i + 1; }
  }
  legs.push({segStart, segEnd: n - 1});
  let prevDir = null, prevLen = 0;
  for (const leg of legs) {
    const p0 = pts[leg.segStart], p1 = pts[leg.segEnd + 1];
    const dir = vnormalize(vsub(p1, p0));
    const len = distance(p0, p1);
    if (prevDir === null) {
      // 最初のレッグ: 対比する入射方向が無いので、自身の延長線上から見返す
      // （光源から近づいてくる光子を正面から見る構図）。
      const offset = Math.min(30, Math.max(5, 0.4 * len));
      leg.aimAt = p1;
      leg.campos = [p1[0] + dir[0]*offset, p1[1] + dir[1]*offset, p1[2] + dir[2]*offset];
    } else {
      const forward = turnForward(prevDir, dir);
      const offset = Math.min(40, Math.max(8, 0.5 * Math.max(len, prevLen)));
      leg.aimAt = p0;
      leg.campos = [p0[0] - forward[0]*offset, p0[1] - forward[1]*offset, p0[2] - forward[2]*offset];
    }
    prevDir = dir; prevLen = len;
  }
  return legs;
}

function legForSeg(legs, seg) {
  for (const leg of legs) if (seg >= leg.segStart && seg <= leg.segEnd) return leg;
  return legs[legs.length - 1];
}

// 区間segが最初に開始した時刻（一人称視点のslerp起点の算出用）。
function buildSegStarts(frames) {
  const starts = [];
  for (const f of frames) if (starts[f.seg] === undefined) starts[f.seg] = f.t0;
  return starts;
}

// 各区間終了時点までの累積反応回数（HUD表示用）
function buildInteractionCounts(traj) {
  const counts = [];
  let c = 0;
  for (const ev of traj.events) { if (INTERACTION_EVENTS.has(ev)) c++; counts.push(c); }
  return counts;
}

let current = null;   // {traj, frames, duration, interactionCounts, idx}
let playState = {clock: 0, playing: false, speed: 1, lastTs: null};
let eMaxGlobal = 0;
for (const t of DATA.trajectories) for (const e of t.energies) if (e > eMaxGlobal) eMaxGlobal = e;
eMaxGlobal = eMaxGlobal || 1;

function energyColor(e) {
  if (!colorEnabled) return '#e2e8f0';
  return `hsl(${240 * (e / eMaxGlobal)},85%,60%)`;
}

function selectPhoton(idx) {
  const traj = DATA.trajectories[idx];
  const {frames, duration} = buildFrames(traj);
  current = {traj, frames, duration, interactionCounts: buildInteractionCounts(traj),
             legs: buildLegs(traj), segStarts: buildSegStarts(frames),
             upAtSegStart: buildUpAtSegStart(traj), idx};
  lastSoundFrame = null;
  playState = {clock: 0, playing: true, speed: playState.speed, lastTs: null};
  document.querySelectorAll('.photon-item').forEach(el =>
    el.classList.toggle('selected', parseInt(el.dataset.idx, 10) === idx));
  document.getElementById('placeholder').style.display = 'none';
  document.getElementById('hud').style.display = 'block';
  document.getElementById('playbar').style.display = 'flex';
  document.getElementById('btnPlay').textContent = '⏸';
  buildRawPanel(traj);
  requestAnimationFrame(tick);
}

// 演出強度設定（スローモー倍率・停止時間）の変更を、選択中の光子の
// 再生タイムラインへ即時反映する。再生位置は経過割合を保って引き継ぐ。
function rebuildCurrentTimeline() {
  if (!current) return;
  const frac = current.duration > 0 ? playState.clock / current.duration : 0;
  const {traj} = current;
  const {frames, duration} = buildFrames(traj);
  current.frames = frames;
  current.duration = duration;
  current.legs = buildLegs(traj);
  current.segStarts = buildSegStarts(frames);
  playState.clock = frac * duration;
  draw();
}

function stateAtClock(clock) {
  const {frames} = current;
  let f = frames[frames.length - 1];
  for (const fr of frames) { if (clock <= fr.t1) { f = fr; break; } }
  const frac = f.t1 > f.t0 ? Math.min(1, Math.max(0, (clock - f.t0) / (f.t1 - f.t0))) : 1;
  const pos = f.pause ? f.from : lerp3(f.from, f.to, frac);
  return {pos, frame: f};
}

function scatterAngleAt(traj, seg) {
  const p0 = traj.points[seg], p1 = traj.points[seg+1], p2 = traj.points[seg+2];
  if (!p2) return null;
  const dIn = [p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]];
  const dOut = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]];
  const nIn = Math.hypot(...dIn), nOut = Math.hypot(...dOut);
  if (nIn === 0 || nOut === 0) return null;
  const cosT = (dIn[0]*dOut[0]+dIn[1]*dOut[1]+dIn[2]*dOut[2]) / (nIn*nOut);
  return Math.acos(Math.min(1, Math.max(-1, cosT))) * 180 / Math.PI;
}

// ---- 生データ表示モード: 演出を介さず、区間ごとの実データをそのまま表形式で見せる ----
function buildRawPanel(traj) {
  const rows = traj.events.map((ev, i) => {
    const len = distance(traj.points[i], traj.points[i+1]);
    const ang = (ev === 'compton' || ev === 'rayleigh') ? scatterAngleAt(traj, i) : null;
    return `<tr data-seg="${i}"><td>${i}</td><td>${traj.materials[i]}</td>` +
           `<td>${traj.energies[i].toFixed(1)}</td><td>${len.toFixed(2)}</td>` +
           `<td>${ev}</td><td>${ang != null ? ang.toFixed(1) + '°' : '—'}</td></tr>`;
  }).join('');
  document.getElementById('rawPanel').innerHTML =
    '<table><thead><tr><th>#</th><th>材料</th><th>E[keV]</th><th>長さ[cm]</th>' +
    '<th>event</th><th>散乱角</th></tr></thead><tbody>' + rows + '</tbody></table>';
}

function updateRawPanelHighlight(seg) {
  const panel = document.getElementById('rawPanel');
  if (panel.style.display !== 'block') return;
  panel.querySelectorAll('tr.current').forEach(tr => tr.classList.remove('current'));
  const row = panel.querySelector(`tr[data-seg="${seg}"]`);
  if (row) row.classList.add('current');
}

function updateHud(state) {
  const {traj} = current;
  const {frame} = state;
  const e = frame.energy, lam = 1239.84 / e;
  const nInter = current.interactionCounts[frame.seg];
  let extra = '';
  if (frame.pause) {
    if (frame.event === 'fluorescence') extra = '<div class="row"><span class="label">事象</span><span>再放出（等方）</span></div>';
    else if (frame.event === 'compton' || frame.event === 'rayleigh') {
      const ang = scatterAngleAt(traj, frame.seg);
      if (ang != null) extra = `<div class="row"><span class="label">散乱角</span><span>${ang.toFixed(1)}°</span></div>`;
    } else if (frame.event === 'photoelectric') {
      extra = '<div class="row"><span class="label">事象</span><span>光電吸収（終了）</span></div>';
    } else if (frame.event === 'escape') {
      extra = '<div class="row"><span class="label">事象</span><span>系外へ脱出（終了）</span></div>';
    }
  }
  document.getElementById('hud').innerHTML =
    `<div class="row"><span class="label">エネルギー</span><span>${e.toFixed(1)} keV</span></div>` +
    `<div class="row"><span class="label">波長</span><span>${lam.toFixed(2)} pm</span></div>` +
    `<div class="row"><span class="label">材料</span><span>${frame.material}</span></div>` +
    `<div class="row"><span class="label">反応回数</span><span>${nInter}</span></div>` + extra;
}

// 一人称視点では終端(吸収/脱出)後に暗転/明転演出のぶんだけ再生時間を延長する。
function effectiveDuration() {
  if (!current) return 0;
  return current.duration + (cameraMode === 'fp' ? TUNE.terminalFadeDuration : 0);
}

function tick(ts) {
  if (!current) return;
  const dur = effectiveDuration();
  if (playState.playing) {
    if (playState.lastTs != null) {
      const dt = (ts - playState.lastTs) / 1000;
      playState.clock = Math.min(dur, playState.clock + dt * playState.speed);
    }
    playState.lastTs = ts;
    if (playState.clock >= dur) {
      playState.playing = false;
      document.getElementById('btnPlay').textContent = '▶';
    }
  } else {
    playState.lastTs = null;
  }
  // シークバーの分母は常にcurrent.duration（カメラモードに依存しない）で固定する。
  // effectiveDuration()を使うと、一人称視点だけ終端フェード秒数ぶん分母が
  // 増えるため、再生中にカメラモードを切り替えた瞬間バー位置が飛んでいた。
  // フェード区間はバー上では「100%で止まったまま裏で再生中」として扱う。
  document.getElementById('seek').value =
    String(Math.min(1000, Math.round(1000 * playState.clock / current.duration)));
  draw();
  requestAnimationFrame(tick);
}

document.getElementById('btnPlay').onclick = () => {
  if (!current) return;
  if (!playState.playing && playState.clock >= effectiveDuration()) playState.clock = 0;
  playState.playing = !playState.playing;
  playState.lastTs = null;
  document.getElementById('btnPlay').textContent = playState.playing ? '⏸' : '▶';
};
document.getElementById('btnRestart').onclick = () => {
  if (!current) return;
  playState.clock = 0; playState.playing = true; playState.lastTs = null;
  document.getElementById('btnPlay').textContent = '⏸';
};
document.getElementById('seek').oninput = e => {
  if (!current) return;
  playState.playing = false; playState.lastTs = null;
  document.getElementById('btnPlay').textContent = '▶';
  playState.clock = current.duration * (parseInt(e.target.value, 10) / 1000);
  draw();
};
document.querySelectorAll('#speedGroup button').forEach(b => b.onclick = () => {
  document.querySelectorAll('#speedGroup button').forEach(x => x.classList.remove('active'));
  b.classList.add('active');
  playState.speed = parseFloat(b.dataset.speed);
});

// 一人称視点の演出オーバーレイ: 大角度転換時の軽いフェード(周辺減光)と、
// 終端(吸収/脱出)後1秒かけての暗転/明転＋メッセージ。物理量ではなく
// 演出であることは画面隅の常設注記(#note)で別途明示している。
function drawFirstPersonOverlay(traj, fpResult, clock, w, h) {
  if (fpResult && fpResult.vignette > 0) {
    const g = ctx.createRadialGradient(w/2, h/2, Math.min(w,h)*0.25, w/2, h/2, Math.min(w,h)*0.7);
    g.addColorStop(0, 'rgba(0,0,0,0)');
    g.addColorStop(1, `rgba(0,0,0,${fpResult.vignette})`);
    ctx.fillStyle = g; ctx.fillRect(0, 0, w, h);
  }
  if (clock > current.duration) {
    const t = Math.min(1, (clock - current.duration) / TUNE.terminalFadeDuration);
    const absorbed = traj.summary.fate === 'absorbed';
    ctx.fillStyle = absorbed ? `rgba(0,0,0,${t})` : `rgba(255,255,255,${t})`;
    ctx.fillRect(0, 0, w, h);
    if (t > 0.5) {
      ctx.globalAlpha = Math.min(1, (t - 0.5) * 2);
      ctx.fillStyle = absorbed ? '#e2e8f0' : '#0f172a';
      ctx.font = 'bold 16px sans-serif';
      ctx.textAlign = 'center';
      const msg = absorbed ? '光電吸収 — 光子はここで消滅した'
        : `系外へ脱出（E = ${traj.energies[traj.energies.length-1].toFixed(1)} keV）`;
      ctx.fillText(msg, w/2, h/2);
      ctx.textAlign = 'left';
      ctx.globalAlpha = 1;
    }
  }
}

function draw() {
  const w = cv.clientWidth, h = cv.clientHeight;

  // activeCamは描画コード全体(drawAxes/drawScene含む)が参照するため、
  // 何を描く前に確定させる。対向リレー視点では、カメラ姿勢はレッグ内で
  // 完全に固定（leg.aimAtという固定点を注視、光子の現在位置は追わない
  // ——buildLegsのコメント参照。追ってしまうと入射側がすぐ視野外に出る）。
  // レッグが変わった瞬間にカメラ位置がカットで切り替わる。
  let state = null, fpResult = null;
  if (current) {
    state = stateAtClock(playState.clock);
    const curSeg = state.frame.seg;
    if (cameraMode === 'relay') {
      const leg = legForSeg(current.legs, curSeg);
      activeCam = buildBasis(leg.campos, leg.aimAt);
    } else if (cameraMode === 'fp') {
      fpResult = computeFirstPersonCam(current.traj, current.segStarts, current.upAtSegStart,
                                        playState.clock, curSeg, state.pos);
      activeCam = fpResult.cam;
    } else {
      activeCam = null;
    }
  } else {
    activeCam = null;
  }

  ctx.clearRect(0, 0, w, h);
  drawAxes();
  drawScene();

  if (current) {
    const {traj} = current;
    const curSeg = state.frame.seg;
    // 通過済みの区間はそのまま描画
    for (let i = 0; i < curSeg; i++) {
      seg(traj.points[i], traj.points[i+1], energyColor(traj.energies[i]), 1.6);
      if (traj.events[i] !== 'boundary')
        markerAt(traj.points[i+1], traj.events[i], energyColor(traj.energies[i]));
    }
    // 現在区間は光子位置まで
    seg(traj.points[curSeg], state.pos, energyColor(traj.energies[curSeg]), 2.2);
    if (state.frame.pause) markerAt(state.pos, state.frame.event, energyColor(state.frame.energy));
    dot(state.pos, 5, '#ffffff');
    dot(state.pos, 3, energyColor(state.frame.energy));
    updateHud(state);
    updateRawPanelHighlight(curSeg);
    if (state.frame.pause && state.frame !== lastSoundFrame) {
      playBlip(state.frame.energy);
      lastSoundFrame = state.frame;
    }

    if (cameraMode === 'fp') drawFirstPersonOverlay(traj, fpResult, playState.clock, w, h);
  }
}

buildSidePanel();
resize();
</script>
</body>
</html>
"""
