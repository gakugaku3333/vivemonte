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
        spectrum = src.get("spectrum")
        if spectrum is not None:
            # 陽なスペクトル指定（単色ビーム等）。kvpベースのSpekPy生成を上書きする
            # ため、kvp自体は不要（sample_spectrumがspectrumを優先する、spectrum.py参照）。
            if not (isinstance(spectrum, list) and spectrum
                    and all(isinstance(s, dict)
                            and isinstance(s.get("energy_keV"), (int, float)) and s["energy_keV"] > 0
                            and isinstance(s.get("weight"), (int, float)) and s["weight"] > 0
                            for s in spectrum)):
                errors.append(SceneError("source.spectrum",
                              "spectrum は [{energy_keV: 正の数値, weight: 正の数値}, ...] の形式で"
                              "指定してください（単色ビームなら1要素、例: "
                              "[{energy_keV: 60, weight: 1.0}]）"))
            # spectrum指定時にkvpも残っていると、実輸送(sample_spectrum)はspectrumを
            # 使うのにpreview表示はkvp側を見てしまう食い違いの温床になる（監査所見）。
            # 曖昧さを許さず、spectrumのみの指定を必須にする。
            if kvp is not None:
                errors.append(SceneError("source",
                              "source.spectrum と source.kvp は併用できません（spectrumが実輸送で"
                              "優先されkvpは無視されるため紛らわしいです）。spectrumのみ指定してください"))
            # spectrumを上書きしても、mAs絶対校正（photon_count_through_field）と
            # ヒール効果（heel_spectra_for_source）はkvpベースのSpekPy計算に固定されて
            # おり、spectrumで指定した実際のエネルギーとは食い違う。相対値[Gy/history]
            # のみのユースケースに限定する。
            if src.get("mas") is not None:
                errors.append(SceneError("source",
                              "source.spectrum と source.mas は併用できません（mAs絶対校正は"
                              "kvpベースのSpekPyフルエンスに固定されており、spectrumで上書きした"
                              "実際のエネルギーと食い違います）。相対値[Gy/history]のみで良ければ"
                              "masを外してください"))
            if src.get("ctdi_vol_mGy") is not None:
                errors.append(SceneError("source",
                              "source.spectrum と source.ctdi_vol_mGy は併用できません（未対応の組み合わせ）"))
            if src.get("heel_effect"):
                errors.append(SceneError("source",
                              "source.spectrum と source.heel_effect は併用できません（ヒール効果は"
                              "kvpベースのSpekPy軸外スペクトルに固定されており、spectrumで上書きした"
                              "実際のエネルギーには対応していません）"))
        elif not isinstance(kvp, (int, float)) or not (20 <= kvp <= 200):
            errors.append(SceneError("source.kvp",
                          f"管電圧 kvp={kvp!r} — 診断領域として 20〜200 kV の数値を指定してください"
                          "（単色ビーム等、kvpを使わない場合は source.spectrum を指定してください）"))
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
            fshape = fld.get("shape", "rect")
            fld["shape"] = fshape
            if fshape == "cone":
                dia = fld.get("diameter_cm")
                if not isinstance(dia, (int, float)) or dia <= 0:
                    errors.append(SceneError("source.field.diameter_cm",
                                  "cone照射野にはSID面での開口直径 diameter_cm（正の数値）が必要です"))
            elif fshape in ("rect", "parallel"):
                size = fld.get("size_cm")
                if not (isinstance(size, (list, tuple)) and len(size) == 2
                        and all(isinstance(x, (int, float)) and x > 0 for x in size)):
                    errors.append(SceneError("source.field.size_cm",
                                  "照射野サイズは正の数値2要素 [幅, 高さ] で指定してください"))
            else:
                errors.append(SceneError("source.field.shape",
                              f"shape={fshape!r} — rect（矩形、既定）/ cone（円錐、立体角一様）/ "
                              "parallel（平行ビーム、非発散）のいずれかです"))
            if fshape == "parallel":
                # 平行ビームは発散を持たないためSIDは無意味（未指定でよい）。
                # 絶対線量校正（mas/heel_effect）はSID基準のSpekPy計算に固定されて
                # おり非対応のため、併用を禁止する。
                if src.get("mas") is not None:
                    errors.append(SceneError("source",
                                  "source.field.shape: parallel と source.mas は併用できません"
                                  "（mAs絶対校正はSID基準のSpekPyフルエンスに固定されており、"
                                  "非発散ビームには対応していません）。相対値[Gy/history]のみで"
                                  "良ければmasを外してください"))
                if src.get("heel_effect"):
                    errors.append(SceneError("source",
                                  "source.field.shape: parallel と source.heel_effect は併用できません"
                                  "（ヒール効果はSID基準の軸外スペクトルに固定されています）"))
            else:
                sid = fld.get("sid_cm")
                if not isinstance(sid, (int, float)) or sid <= 0:
                    errors.append(SceneError("source.field.sid_cm", "SID（焦点-照射野定義面距離）は正の数値です"))
        filt = src.get("filtration_mm_al", 2.5)
        if not isinstance(filt, (int, float)) or filt < 0:
            errors.append(SceneError("source.filtration_mm_al", "総濾過は0以上の数値（mmAl）です"))
        elif filt < 1.5 and isinstance(kvp, (int, float)) and kvp >= 70:
            warnings.append(SceneError("source.filtration_mm_al",
                          f"総濾過 {filt} mmAl は診断装置の法令要件（一般に2.5 mmAl以上）より薄い可能性があります"))

        angle = src.get("anode_angle_deg", 12.0)
        if not isinstance(angle, (int, float)) or not (5.0 <= angle <= 20.0):
            warnings.append(SceneError("source.anode_angle_deg",
                          f"陽極角 {angle!r} — 一般的な診断用X線管は5〜20度です（既定12度を使用、SpekPyスペクトル計算に使用）"))
            angle = 12.0
        src["anode_angle_deg"] = float(angle)

        mas = src.get("mas")
        if mas is not None:
            if not isinstance(mas, (int, float)) or mas <= 0:
                errors.append(SceneError("source.mas", "mas（管電流時間積）は正の数値です"))
            else:
                src["mas"] = float(mas)

        ctdi = src.get("ctdi_vol_mGy")
        if ctdi is not None:
            if not isinstance(ctdi, (int, float)) or ctdi <= 0:
                errors.append(SceneError("source.ctdi_vol_mGy",
                              "ctdi_vol_mGy（コンソール表示のCTDIvol）は正の数値です"))
            else:
                src["ctdi_vol_mGy"] = float(ctdi)
            if mas is not None:
                errors.append(SceneError("source",
                              "mas と ctdi_vol_mGy は同時に指定できません（絶対線量校正の基準が"
                              "曖昧になります）。CTなら ctdi_vol_mGy、一般撮影なら mas を使ってください"))
            if src.get("rotation") is None:
                errors.append(SceneError("source.ctdi_vol_mGy",
                              "CTDIvol校正には source.rotation（ガントリー回転）が必要です"))
            ph = src.get("ctdi_phantom", "body")
            src["ctdi_phantom"] = ph
            if ph not in ("body", "head"):
                errors.append(SceneError("source.ctdi_phantom",
                              f"ctdi_phantom={ph!r} — body（Ø32cm）/ head（Ø16cm）のいずれかです"))

        if src.get("heel_effect"):
            ad = _vec3(src.get("anode_direction"), "source.anode_direction", errors,
                       "陽極方向（heel_effect使用時は必須）")
            if ad is not None and isinstance(src.get("direction"), list):
                dv = src["direction"]  # 上で正規化済み
                perp = [ad[k] - (ad[0] * dv[0] + ad[1] * dv[1] + ad[2] * dv[2]) * dv[k]
                        for k in range(3)]
                if math.sqrt(sum(x * x for x in perp)) < 1e-6:
                    errors.append(SceneError("source.anode_direction",
                                  "陽極方向がビーム中心軸と平行です。陽極-陰極軸は"
                                  "中心軸に直交する成分を持つ必要があります"))

        rot = src.get("rotation")
        if rot is not None:
            if not isinstance(rot, dict):
                errors.append(SceneError("source.rotation", "rotation はマッピングで指定してください"))
            else:
                _vec3(rot.get("isocenter"), "source.rotation.isocenter", errors, "回転中心座標")
                axis = rot.get("axis", "z")
                rot["axis"] = axis
                if axis not in VALID_AXES:
                    errors.append(SceneError("source.rotation.axis", f"axis={axis!r} — x/y/zのいずれかです"))
                n_angles = rot.get("n_angles")
                if n_angles is not None and not (isinstance(n_angles, int) and n_angles >= 2):
                    errors.append(SceneError("source.rotation.n_angles",
                                  "n_angles（離散ガントリー角度数）は2以上の整数です。"
                                  "省略時は連続一様（CTの連続曝射に対応、通常はこちらを推奨）"))
                scan = rot.get("scan_length_cm")
                if scan is not None:
                    if not isinstance(scan, (int, float)) or scan <= 0:
                        errors.append(SceneError("source.rotation.scan_length_cm",
                                      "scan_length_cm（ヘリカルスキャン範囲）は正の数値です"))
                    else:
                        rot["scan_length_cm"] = float(scan)

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
    """照射野定義面（SID位置。parallelはposition自体の面）における照射野の外周点を返す。

    rect/parallel照射野は4隅、cone照射野は開口円周上の16点（描画用の多角形近似）。
    """
    pos = src["position"]
    d = src["direction"]
    fld = src["field"]
    shape = fld.get("shape", "rect")
    # 中心軸に直交する基底ベクトル（uを水平寄り、vをその直交に取る）
    if abs(d[2]) < 0.999:
        u = [-d[1], d[0], 0.0]
    else:
        u = [1.0, 0.0, 0.0]
    n = math.sqrt(sum(x * x for x in u))
    u = [x / n for x in u]
    v = [d[1] * u[2] - d[2] * u[1], d[2] * u[0] - d[0] * u[2], d[0] * u[1] - d[1] * u[0]]
    if shape == "parallel":
        # 非発散ビーム: 定義面はposition自体（SIDを持たない）
        w, h = fld["size_cm"]
        return [[pos[k] + su * w / 2 * u[k] + sv * h / 2 * v[k] for k in range(3)]
                for su, sv in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    ctr = [pos[k] + d[k] * fld["sid_cm"] for k in range(3)]
    if shape == "cone":
        r = fld["diameter_cm"] / 2.0
        pts = []
        for i in range(16):
            th = 2 * math.pi * i / 16
            pts.append([ctr[k] + r * math.cos(th) * u[k] + r * math.sin(th) * v[k]
                        for k in range(3)])
        return pts
    w, h = fld["size_cm"]
    corners = []
    for su, sv in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        corners.append([ctr[k] + su * w / 2 * u[k] + sv * h / 2 * v[k] for k in range(3)])
    return corners
