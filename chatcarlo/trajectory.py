"""光子軌跡の記録と可視化用データ整形（`chatcarlo trace` の裏方）。

輸送カーネル（chatcarlo/transport.py）に TrajectoryRecorder を渡すと、
ループ1周ごとの飛行区間が記録される。乱数を一切消費しないため、
recorderの有無は同一seedでの輸送結果に影響しない。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class TrajectoryRecorder:
    """軌跡記録（小history可視化用）。ループ1周ごとに飛行区間を追記する。

    starts/ends/energies/events/materials/photon_ids はそれぞれ「1反復ぶんの配列」の
    リストとして貯め、trajectories_to_json() で光子ごとのポリラインにまとめる。
    event は区間の終端で起きたことを表す文字列:
      "boundary"（材料境界を通過して継続）, "photoelectric", "compton",
      "rayleigh", "fluorescence"（K殻蛍光X線放出、光電の亜種）, "escape"
    material は区間**始点**での材料名（geometry.material_at(o)、区間内で材料は
    一定なので終点でも同じ）。
    """
    starts: list = field(default_factory=list)
    ends: list = field(default_factory=list)
    energies: list = field(default_factory=list)
    events: list = field(default_factory=list)
    materials: list = field(default_factory=list)
    photon_ids: list = field(default_factory=list)

    def record(self, photon_id: np.ndarray, start: np.ndarray, end: np.ndarray,
               energy_keV: np.ndarray, event: np.ndarray,
               material: np.ndarray | None = None) -> None:
        self.photon_ids.append(np.asarray(photon_id))
        self.starts.append(np.asarray(start))
        self.ends.append(np.asarray(end))
        self.energies.append(np.asarray(energy_keV, dtype=float))
        self.events.append(np.asarray(event, dtype=object))
        if material is None:
            material = np.full(len(photon_id), "", dtype=object)
        self.materials.append(np.asarray(material, dtype=object))


def trajectories_to_json(recorder: TrajectoryRecorder) -> list[dict]:
    """TrajectoryRecorderの飛行区間データを光子ごとのポリラインにまとめる。

    区間はrecorderへの追記順（=輸送ループの反復順）であり、同一photon_idの
    区間は反復ごとに高々1つしか記録されないため、そのまま連結すれば
    時系列順のポリラインになる。
    """
    if not recorder.photon_ids:
        return []
    photon_ids = np.concatenate(recorder.photon_ids)
    starts = np.concatenate(recorder.starts)
    ends = np.concatenate(recorder.ends)
    energies = np.concatenate(recorder.energies)
    events = np.concatenate(recorder.events)
    materials = np.concatenate(recorder.materials)

    by_photon: dict[int, dict] = {}
    order: list[int] = []
    for i in range(len(photon_ids)):
        pid = int(photon_ids[i])
        traj = by_photon.get(pid)
        if traj is None:
            traj = {"points": [starts[i].tolist()], "energies": [], "events": [], "materials": []}
            by_photon[pid] = traj
            order.append(pid)
        traj["points"].append(ends[i].tolist())
        traj["energies"].append(float(energies[i]))
        traj["events"].append(str(events[i]))
        traj["materials"].append(str(materials[i]))

    trajectories = [by_photon[pid] for pid in order]
    for traj in trajectories:
        traj["summary"] = classify_trajectory(traj)
    return trajectories


def classify_trajectory(traj: dict) -> dict:
    """光子1本の転帰分類（trajectories_to_json の1要素を受け取る）。

    events は区間の終端理由の列（"boundary"/"photoelectric"/"compton"/
    "rayleigh"/"fluorescence"/"escape"）。transport.pyのループ終了条件は
    「光電吸収（蛍光放出なし）」か「脱出」の2つだけなので、最終eventは
    必ず "photoelectric" か "escape" のいずれかになる。

    戻り値:
      n_interactions: photoelectric/compton/rayleigh/fluorescence の件数
        （fluorescenceは光電吸収イベントの亜種＝再放出であり、コンプトン/
        レイリーの散乱とは運動学が異なるが、ここでは「反応が起きた区間」
        として同じく1件に数える）
      fate: "absorbed"（最終eventがphotoelectric、蛍光なしで局所吸収）
            | "escaped"（最終eventがescape）
      backscattered: 初期方向と最終区間の方向の内積が負ならTrue
      category: "direct"（無散乱で脱出）| "scattered_escape"（散乱後に脱出）
                | "backscatter"（脱出かつbackscattered）| "absorbed"
    """
    events = traj["events"]
    points = traj["points"]
    interaction_events = {"photoelectric", "compton", "rayleigh", "fluorescence"}
    n_interactions = sum(1 for ev in events if ev in interaction_events)
    fate = "absorbed" if events[-1] == "photoelectric" else "escaped"

    p0, p1 = np.asarray(points[0]), np.asarray(points[1])
    pm1, pm2 = np.asarray(points[-1]), np.asarray(points[-2])
    d_init = p1 - p0
    d_final = pm1 - pm2
    n_init = np.linalg.norm(d_init)
    n_final = np.linalg.norm(d_final)
    backscattered = bool(
        n_init > 0 and n_final > 0
        and np.dot(d_init, d_final) / (n_init * n_final) < 0
    )

    if fate == "absorbed":
        category = "absorbed"
    elif backscattered:
        category = "backscatter"
    elif n_interactions == 0:
        category = "direct"
    else:
        category = "scattered_escape"

    return {"n_interactions": n_interactions, "fate": fate,
            "backscattered": backscattered, "category": category}
