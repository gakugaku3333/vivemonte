"""ジオメトリープレビュー — scene.yaml から自己完結HTML 3Dビューアを生成する。

外部ライブラリ・CDNに依存しない（vanilla JS + canvas）。
AIが組んだ体系をユーザーが目視確認するための、承認ゲートの中核。
"""
from __future__ import annotations

import json
import math

from .scene import Scene, field_corners

# 材料ごとの表示色（見分けやすさ優先）
MATERIAL_COLORS = {
    "water": "#3b82f6", "air": "#94a3b8", "soft_tissue": "#f59e0b",
    "bone": "#e2e8f0", "lung": "#f9a8d4", "pmma": "#22d3ee",
    "concrete": "#a8a29e", "lead": "#8b5cf6", "aluminum": "#67e8f9",
    "copper": "#fb923c", "iron": "#f87171", "lead_glass": "#c084fc",
}
_FALLBACK = ["#10b981", "#eab308", "#ec4899", "#14b8a6", "#f97316"]


def _box_mesh(c, s):
    hx, hy, hz = s[0] / 2, s[1] / 2, s[2] / 2
    v = [[c[0] + sx * hx, c[1] + sy * hy, c[2] + sz * hz]
         for sz in (-1, 1) for sy in (-1, 1) for sx in (-1, 1)]
    edges_i = [(0, 1), (1, 3), (3, 2), (2, 0), (4, 5), (5, 7), (7, 6), (6, 4),
               (0, 4), (1, 5), (2, 6), (3, 7)]
    faces_i = [(0, 1, 3, 2), (4, 5, 7, 6), (0, 1, 5, 4), (2, 3, 7, 6),
               (0, 2, 6, 4), (1, 3, 7, 5)]
    return ([[v[a], v[b]] for a, b in edges_i],
            [[v[i] for i in f] for f in faces_i])


def _cylinder_mesh(c, r, h, axis, n=24):
    ax = {"x": 0, "y": 1, "z": 2}[axis]
    u, w = [i for i in range(3) if i != ax]
    rings = []
    for end in (-1, 1):
        ring = []
        for k in range(n):
            th = 2 * math.pi * k / n
            p = list(c)
            p[ax] += end * h / 2
            p[u] += r * math.cos(th)
            p[w] += r * math.sin(th)
            ring.append(p)
        rings.append(ring)
    bot, top = rings
    edges = []
    for ring in rings:
        edges += [[ring[k], ring[(k + 1) % n]] for k in range(n)]
    edges += [[bot[k], top[k]] for k in range(0, n, 3)]
    faces = [[bot[k], bot[(k + 1) % n], top[(k + 1) % n], top[k]] for k in range(n)]
    faces += [bot, top]
    return edges, faces


def _sphere_mesh(c, r, n=24):
    edges = []
    for plane in ((0, 1), (0, 2), (1, 2)):
        ring = []
        for k in range(n):
            th = 2 * math.pi * k / n
            p = list(c)
            p[plane[0]] += r * math.cos(th)
            p[plane[1]] += r * math.sin(th)
            ring.append(p)
        edges += [[ring[k], ring[(k + 1) % n]] for k in range(n)]
    return edges, []


def scene_to_json(scene: Scene, trajectories: list[dict] | None = None) -> dict:
    raw = scene.raw
    objects = []
    mats_seen = {}
    for i, g in enumerate(raw.get("geometry", [])):
        mat = str(g.get("material", "?"))
        if mat not in mats_seen:
            mats_seen[mat] = MATERIAL_COLORS.get(mat.lower(),
                                                 _FALLBACK[len(mats_seen) % len(_FALLBACK)])
        color = g.get("color", mats_seen[mat])
        shape, c = g["shape"], g["center"]
        if shape == "box":
            edges, faces = _box_mesh(c, g["size_cm"])
        elif shape == "cylinder":
            edges, faces = _cylinder_mesh(c, g["radius_cm"], g["height_cm"], g.get("axis", "z"))
        else:
            edges, faces = _sphere_mesh(c, g["radius_cm"])
        objects.append({"name": g["name"], "material": mat, "color": color,
                        "edges": edges, "faces": faces, "label_at": c})

    src = raw["source"]
    corners = field_corners(src)
    pos = src["position"]
    beam = {"source": pos, "corners": corners, "kvp": src["kvp"],
            "field_size": src["field"]["size_cm"], "sid": src["field"]["sid_cm"]}

    # シーン全体のバウンディング
    pts = [pos] + corners
    for o in objects:
        for e in o["edges"]:
            pts += e
    lo = [min(p[k] for p in pts) for k in range(3)]
    hi = [max(p[k] for p in pts) for k in range(3)]
    center = [(lo[k] + hi[k]) / 2 for k in range(3)]
    radius = max(math.dist(lo, hi) / 2, 1.0)

    return {"objects": objects, "beam": beam, "center": center, "radius": radius,
            "warnings": [str(w) for w in scene.warnings], "trajectories": trajectories or []}


def render_html(scene: Scene, title: str = "viveMonte geometry preview",
                 trajectories: list[dict] | None = None) -> str:
    data = json.dumps(scene_to_json(scene, trajectories=trajectories), ensure_ascii=False)
    return _TEMPLATE.replace("__TITLE__", title).replace("__DATA__", data)


def write_html(scene: Scene, out_path: str, title: str = "viveMonte geometry preview",
                trajectories: list[dict] | None = None) -> str:
    html = render_html(scene, title, trajectories=trajectories)
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
  header { padding:10px 16px; display:flex; gap:12px; align-items:center; flex-wrap:wrap; }
  header h1 { font-size:15px; font-weight:600; margin-right:auto; }
  button { background:var(--panel); color:var(--text); border:1px solid #47556933;
           border-radius:6px; padding:5px 12px; font-size:12px; cursor:pointer; }
  button:hover { border-color:var(--accent); }
  label { font-size:12px; color:var(--dim); display:flex; align-items:center; gap:4px; cursor:pointer; }
  #wrap { flex:1; position:relative; }
  canvas { position:absolute; inset:0; width:100%; height:100%; cursor:grab; }
  #legend { position:absolute; top:10px; left:12px; background:var(--panel); border-radius:8px;
            padding:10px 14px; font-size:12px; line-height:1.9; opacity:.94; max-width:260px; }
  #legend .sw { display:inline-block; width:11px; height:11px; border-radius:2px; margin-right:7px; }
  #info { position:absolute; bottom:10px; left:12px; font-size:11px; color:var(--dim); }
  #warn { position:absolute; top:10px; right:12px; max-width:320px; background:#7c2d12; color:#fed7aa;
          border-radius:8px; padding:8px 12px; font-size:12px; line-height:1.5; }
</style>
</head>
<body>
<header>
  <h1>__TITLE__</h1>
  <button data-view="iso">等角</button>
  <button data-view="front">正面 (−Y)</button>
  <button data-view="side">側面 (+X)</button>
  <button data-view="top">上面 (+Z)</button>
  <label><input type="checkbox" id="ckFaces" checked> 面塗り</label>
  <label><input type="checkbox" id="ckBeam" checked> ビーム</label>
  <label><input type="checkbox" id="ckLabels" checked> ラベル</label>
  <label id="lblTraj" style="display:none"><input type="checkbox" id="ckTraj" checked> 軌跡</label>
</header>
<div id="wrap">
  <canvas id="cv"></canvas>
  <div id="legend"></div>
  <div id="info">ドラッグ: 回転　/　Shift+ドラッグ: 平行移動　/　ホイール: ズーム</div>
  <div id="warn" style="display:none"></div>
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
['ckFaces','ckBeam','ckLabels','ckTraj'].forEach(id => document.getElementById(id).onchange = draw);
if (DATA.trajectories.length) document.getElementById('lblTraj').style.display = 'flex';

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

// 世界座標(z上向き) → カメラ → 画面。カメラは中心を注視して周回。
function project(p) {
  const c = DATA.center;
  let x = p[0]-c[0], y = p[1]-c[1], z = p[2]-c[2];
  const cy = Math.cos(yaw), sy = Math.sin(yaw);
  let x1 = x*cy - y*sy, y1 = x*sy + y*cy;
  const cp = Math.cos(pitch), sp = Math.sin(pitch);
  let y2 = y1*cp - z*sp, z2 = y1*sp + z*cp;      // y2: 奥行き, z2: 上
  const depth = y2 + dist;
  if (depth < 1e-3) return null;
  const w = cv.clientWidth, h = cv.clientHeight;
  const f = 1.2 * Math.min(w, h) / (2 * Math.tan(0.35)) / dist;
  return [w/2 + x1 * f * dist/depth + panX, h/2 - z2 * f * dist/depth + panY, depth];
}

function poly(pts, fill, alpha) {
  const pr = pts.map(project); if (pr.some(p => !p)) return;
  ctx.beginPath(); ctx.moveTo(pr[0][0], pr[0][1]);
  for (let i = 1; i < pr.length; i++) ctx.lineTo(pr[i][0], pr[i][1]);
  ctx.closePath(); ctx.globalAlpha = alpha; ctx.fillStyle = fill; ctx.fill(); ctx.globalAlpha = 1;
}
function seg(a, b, color, width, dash) {
  const pa = project(a), pb = project(b); if (!pa || !pb) return;
  ctx.beginPath(); ctx.moveTo(pa[0], pa[1]); ctx.lineTo(pb[0], pb[1]);
  ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash || []); ctx.stroke(); ctx.setLineDash([]);
}

function drawAxes() {
  const c = DATA.center, L = DATA.radius * 0.55;
  const o = [c[0], c[1], c[2]];
  const ax = [[[L,0,0],'#ef4444','X'], [[0,L,0],'#22c55e','Y'], [[0,0,L],'#3b82f6','Z']];
  for (const [d, col, name] of ax) {
    const e = [o[0]+d[0], o[1]+d[1], o[2]+d[2]];
    seg(o, e, col, 1.5);
    const p = project(e);
    if (p) { ctx.fillStyle = col; ctx.font = '11px sans-serif'; ctx.fillText(name, p[0]+4, p[1]-4); }
  }
}

function draw() {
  const w = cv.clientWidth, h = cv.clientHeight;
  ctx.clearRect(0, 0, w, h);
  const showFaces = ckFaces.checked, showBeam = ckBeam.checked, showLabels = ckLabels.checked;

  drawAxes();

  if (showFaces) {
    // 面を奥から手前へ（painter's algorithm・重心深度）
    const all = [];
    for (const o of DATA.objects)
      for (const f of o.faces) {
        let d = 0, n = 0;
        for (const p of f) { const pr = project(p); if (!pr) { d = -1; break; } d += pr[2]; n++; }
        if (d > 0) all.push([d / n, f, o.color]);
      }
    all.sort((a, b) => b[0] - a[0]);
    for (const [, f, col] of all) poly(f, col, 0.13);
  }
  for (const o of DATA.objects)
    for (const e of o.edges) seg(e[0], e[1], o.color, 1.1);

  if (showBeam) {
    const b = DATA.beam;
    for (const c of b.corners) seg(b.source, c, '#fbbf24', 1.2, [5, 4]);
    for (let i = 0; i < 4; i++) seg(b.corners[i], b.corners[(i+1)%4], '#fbbf24', 2);
    poly(b.corners, '#fbbf24', 0.18);
    const ps = project(b.source);
    if (ps) {
      ctx.fillStyle = '#fbbf24'; ctx.beginPath(); ctx.arc(ps[0], ps[1], 5, 0, 7); ctx.fill();
      ctx.font = 'bold 12px sans-serif';
      ctx.fillText(`焦点 ${b.kvp} kV`, ps[0]+9, ps[1]+4);
    }
  }
  if (ckTraj.checked && DATA.trajectories.length) {
    let eMax = 0;
    for (const t of DATA.trajectories) for (const en of t.energies) if (en > eMax) eMax = en;
    eMax = eMax || 1;
    for (const t of DATA.trajectories) {
      const pts = t.points;
      for (let i = 0; i < t.energies.length; i++) {
        const hue = 240 * (t.energies[i] / eMax);
        const col = `hsl(${hue},85%,60%)`;
        seg(pts[i], pts[i + 1], col, 1.4);
        const p = project(pts[i + 1]);
        if (!p) continue;
        const ev = t.events[i];
        if (ev === 'photoelectric') {
          ctx.fillStyle = col; ctx.beginPath(); ctx.arc(p[0], p[1], 3, 0, 7); ctx.fill();
        } else if (ev === 'compton') {
          ctx.strokeStyle = col; ctx.lineWidth = 1;
          ctx.beginPath(); ctx.arc(p[0], p[1], 3, 0, 7); ctx.stroke();
        } else if (ev === 'rayleigh') {
          ctx.strokeStyle = col; ctx.lineWidth = 1;
          ctx.beginPath();
          ctx.moveTo(p[0], p[1] - 4); ctx.lineTo(p[0] + 4, p[1]);
          ctx.lineTo(p[0], p[1] + 4); ctx.lineTo(p[0] - 4, p[1]);
          ctx.closePath(); ctx.stroke();
        } else if (ev === 'escape') {
          ctx.strokeStyle = col; ctx.lineWidth = 1.3;
          ctx.beginPath();
          ctx.moveTo(p[0] - 3, p[1] - 3); ctx.lineTo(p[0] + 3, p[1] + 3);
          ctx.moveTo(p[0] - 3, p[1] + 3); ctx.lineTo(p[0] + 3, p[1] - 3);
          ctx.stroke();
        }
      }
    }
  }
  if (showLabels) {
    ctx.font = '12px sans-serif';
    for (const o of DATA.objects) {
      const p = project(o.label_at); if (!p) continue;
      ctx.fillStyle = o.color; ctx.fillText(o.name, p[0]+5, p[1]);
    }
  }
}

// 凡例と警告
{
  const lg = document.getElementById('legend');
  const mats = {};
  for (const o of DATA.objects) mats[o.material] = o.color;
  lg.innerHTML = '<b style="font-size:12px">材料</b><br>' +
    Object.entries(mats).map(([m, c]) =>
      `<span class="sw" style="background:${c}"></span>${m}`).join('<br>') +
    `<br><span class="sw" style="background:#fbbf24"></span>X線ビーム（照射野 ${DATA.beam.field_size[0]}×${DATA.beam.field_size[1]} cm @ SID ${DATA.beam.sid} cm）`;
  if (DATA.trajectories.length) {
    lg.innerHTML += '<br><br><b style="font-size:12px">軌跡（' + DATA.trajectories.length + '光子）</b><br>' +
      '● 光電吸収　○ コンプトン　◇ レイリー　× 脱出<br>' +
      '色: 青=高エネルギー → 赤=低エネルギー';
  }
  if (DATA.warnings.length) {
    const wd = document.getElementById('warn');
    wd.style.display = 'block';
    wd.innerHTML = '<b>検証警告</b><br>' + DATA.warnings.join('<br>');
  }
}
resize();
</script>
</body>
</html>
"""
