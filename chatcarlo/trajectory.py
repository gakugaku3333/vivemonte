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

    starts/ends/energies/events/photon_ids はそれぞれ「1反復ぶんの配列」の
    リストとして貯め、trajectories_to_json() で光子ごとのポリラインにまとめる。
    event は区間の終端で起きたことを表す文字列:
      "boundary"（材料境界を通過して継続）, "photoelectric", "compton",
      "rayleigh", "fluorescence"（K殻蛍光X線放出、光電の亜種）, "escape"
    """
    starts: list = field(default_factory=list)
    ends: list = field(default_factory=list)
    energies: list = field(default_factory=list)
    events: list = field(default_factory=list)
    photon_ids: list = field(default_factory=list)

    def record(self, photon_id: np.ndarray, start: np.ndarray, end: np.ndarray,
               energy_keV: np.ndarray, event: np.ndarray) -> None:
        self.photon_ids.append(np.asarray(photon_id))
        self.starts.append(np.asarray(start))
        self.ends.append(np.asarray(end))
        self.energies.append(np.asarray(energy_keV, dtype=float))
        self.events.append(np.asarray(event, dtype=object))


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

    by_photon: dict[int, dict] = {}
    order: list[int] = []
    for i in range(len(photon_ids)):
        pid = int(photon_ids[i])
        traj = by_photon.get(pid)
        if traj is None:
            traj = {"points": [starts[i].tolist()], "energies": [], "events": []}
            by_photon[pid] = traj
            order.append(pid)
        traj["points"].append(ends[i].tolist())
        traj["energies"].append(float(energies[i]))
        traj["events"].append(str(events[i]))

    return [by_photon[pid] for pid in order]
