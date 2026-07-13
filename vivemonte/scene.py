"""scene.yaml の読み込みと検証。

AIが生成したシーン記述を機械検証し、明確なエラーメッセージで
自己修正ループを回せるようにするのがこのモジュールの役割。
座標系: cm単位、z軸が鉛直上向き、床が z=0。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import yaml

VALID_SHAPES = {"box", "cylinder", "sphere"}
VALID_AXES = {"x", "y", "z"}


@dataclass
class SceneError:
    path: str      # 例: "geometry[2].size_cm"
    message: str

    def __str__(self) -> str:
        return f"{self.path}: {self.message}"


@dataclass
class Scene:
    raw: dict
    errors: list[SceneError] = field(default_factory=list)
    warnings: list[SceneError] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def _vec3(v, path, errors, name="ベクトル"):
    if not (isinstance(v, (list, tuple)) and len(v) == 3
            and all(isinstance(x, (int, float)) for x in v)):
        errors.append(SceneError(path, f"{name}は数値3要素のリストで指定してください（例: [0, 0, 100]）"))
        return None
    return [float(x) for x in v]


def load_scene(path: str) -> Scene:
    with open(path, encoding="utf-8") as f:
        raw = yaml.safe_load(f)
    return validate_scene(raw)


def validate_scene(raw: dict) -> Scene:
    errors: list[SceneError] = []
    warnings: list[SceneError] = []
    scene = Scene(raw=raw, errors=errors, warnings=warnings)

    if not isinstance(raw, dict):
        errors.append(SceneError("(root)", "YAMLのトップレベルはマッピングである必要があります"))
        return scene

    # ---- source ----
    src = raw.get("source")
    if src is None:
        errors.append(SceneError("source", "source セクションがありません"))
    else:
        kvp = src.get("kvp")
        if not isinstance(kvp, (int, float)) or not (20 <= kvp <= 200):
            errors.append(SceneError("source.kvp",
                          f"管電圧 kvp={kvp!r} — 診断領域として 20〜200 kV の数値を指定してください"))
        pos = _vec3(src.get("position"), "source.position", errors, "焦点位置")
        dirv = _vec3(src.get("direction"), "source.direction", errors, "中心軸方向")
        if dirv is not None:
            n = math.sqrt(sum(x * x for x in dirv))
            if n < 1e-9:
                errors.append(SceneError("source.direction", "方向ベクトルがゼロです"))
            else:
                src["direction"] = [x / n for x in dirv]
        fld = src.get("field")
        if fld is None:
            errors.append(SceneError("source.field", "照射野 field がありません"))
        else:
            size = fld.get("size_cm")
            if not (isinstance(size, (list, tuple)) and len(size) == 2
                    and all(isinstance(x, (int, float)) and x > 0 for x in size)):
                errors.append(SceneError("source.field.size_cm",
                              "照射野サイズは正の数値2要素 [幅, 高さ] で指定してください"))
            sid = fld.get("sid_cm")
            if not isinstance(sid, (int, float)) or sid <= 0:
                errors.append(SceneError("source.field.sid_cm", "SID（焦点-照射野定義面距離）は正の数値です"))
        filt = src.get("filtration_mm_al", 2.5)
        if not isinstance(filt, (int, float)) or filt < 0:
            errors.append(SceneError("source.filtration_mm_al", "総濾過は0以上の数値（mmAl）です"))
        elif filt < 1.5 and isinstance(kvp, (int, float)) and kvp >= 70:
            warnings.append(SceneError("source.filtration_mm_al",
                          f"総濾過 {filt} mmAl は診断装置の法令要件（一般に2.5 mmAl以上）より薄い可能性があります"))

    # ---- geometry ----
    geoms = raw.get("geometry")
    if not isinstance(geoms, list) or not geoms:
        errors.append(SceneError("geometry", "geometry には1つ以上の物体をリストで指定してください"))
        geoms = []

    names = set()
    for i, g in enumerate(geoms):
        p = f"geometry[{i}]"
        if not isinstance(g, dict):
            errors.append(SceneError(p, "各物体はマッピングで指定してください"))
            continue
        name = g.get("name", f"object_{i}")
        g["name"] = name
        if name in names:
            errors.append(SceneError(f"{p}.name", f"物体名 '{name}' が重複しています"))
        names.add(name)

        shape = g.get("shape")
        if shape not in VALID_SHAPES:
            errors.append(SceneError(f"{p}.shape",
                          f"shape={shape!r} — {sorted(VALID_SHAPES)} のいずれかを指定してください"))
            continue
        if not g.get("material"):
            errors.append(SceneError(f"{p}.material", "material がありません"))

        _vec3(g.get("center"), f"{p}.center", errors, "中心座標")

        if shape == "box":
            size = g.get("size_cm")
            if not (isinstance(size, (list, tuple)) and len(size) == 3
                    and all(isinstance(x, (int, float)) and x > 0 for x in size)):
                errors.append(SceneError(f"{p}.size_cm", "boxは正の数値3要素 [x, y, z] で指定してください"))
        elif shape == "cylinder":
            if not (isinstance(g.get("radius_cm"), (int, float)) and g["radius_cm"] > 0):
                errors.append(SceneError(f"{p}.radius_cm", "cylinderには正の radius_cm が必要です"))
            if not (isinstance(g.get("height_cm"), (int, float)) and g["height_cm"] > 0):
                errors.append(SceneError(f"{p}.height_cm", "cylinderには正の height_cm が必要です"))
            axis = g.get("axis", "z")
            g["axis"] = axis
            if axis not in VALID_AXES:
                errors.append(SceneError(f"{p}.axis", f"axis={axis!r} — x/y/z のいずれかです"))
        elif shape == "sphere":
            if not (isinstance(g.get("radius_cm"), (int, float)) and g["radius_cm"] > 0):
                errors.append(SceneError(f"{p}.radius_cm", "sphereには正の radius_cm が必要です"))

    # ---- 物理サニティチェック（警告） ----
    if scene.ok and src:
        pos = src["position"]
        for g in geoms:
            if g.get("shape") == "box" and "size_cm" in g and "center" in g:
                c, s = g["center"], g["size_cm"]
                inside = all(abs(pos[k] - c[k]) < s[k] / 2 for k in range(3))
                if inside:
                    warnings.append(SceneError("source.position",
                                  f"X線焦点が物体 '{g['name']}' の内部にあります。意図した配置か確認してください"))
    return scene


def field_corners(src: dict) -> list[list[float]]:
    """照射野定義面（SID位置）における照射野4隅の座標を返す。"""
    pos = src["position"]
    d = src["direction"]
    w, h = src["field"]["size_cm"]
    sid = src["field"]["sid_cm"]
    # 中心軸に直交する基底ベクトル（uを水平寄り、vをその直交に取る）
    if abs(d[2]) < 0.999:
        u = [-d[1], d[0], 0.0]
    else:
        u = [1.0, 0.0, 0.0]
    n = math.sqrt(sum(x * x for x in u))
    u = [x / n for x in u]
    v = [d[1] * u[2] - d[2] * u[1], d[2] * u[0] - d[0] * u[2], d[0] * u[1] - d[1] * u[0]]
    ctr = [pos[k] + d[k] * sid for k in range(3)]
    corners = []
    for su, sv in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        corners.append([ctr[k] + su * w / 2 * u[k] + sv * h / 2 * v[k] for k in range(3)])
    return corners
