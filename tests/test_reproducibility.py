"""同一seedの再現性 — プロセス（文字列ハッシュシード）が違っても結果が一致すること。

set()の反復順は文字列ハッシュのランダム化（PYTHONHASHSEED）でプロセスごとに
変わるため、材料グループの順に依存して乱数を消費する処理（レイリー元素抽選
など）が未ソートだと、同一seedでも実行ごとに結果が変わってしまう。
materials.material_groups() のソートがその対策（過去に実際に起きた回帰）。
"""
from __future__ import annotations

import os
import subprocess
import sys

import numpy as np

from vivemonte.materials import material_groups

_SCRIPT = """
import numpy as np
from vivemonte.geometry import Geometry
from vivemonte.source import sample_source_photons
from vivemonte.transport import transport_photons

src = {"kvp": 80, "position": [0, -30, 0], "direction": [0, 1, 0],
       "field": {"size_cm": [10, 10], "sid_cm": 30}}
geometry = Geometry([
    {"name": "slab", "shape": "box", "material": "water",
     "center": [0, 0, 0], "size_cm": [10, 4, 10]},
    {"name": "plate", "shape": "box", "material": "aluminum",
     "center": [0, 5, 0], "size_cm": [10, 1, 10]},
])
rng = np.random.default_rng(7)
pos, dirv, e = sample_source_photons(src, 3000, rng)
r = transport_photons(pos, dirv, e, geometry, rng)
print(sorted((k, round(v, 6)) for k, v in r.energy_deposited.items()))
print(int(r.absorbed.sum()), int(r.escaped.sum()), int(r.n_scatter.sum()))
"""


def test_material_groups_order_is_sorted():
    names = np.array(["water", "aluminum", "water", "lead", "aluminum"])
    got = [name for name, _ in material_groups(names)]
    assert got == ["aluminum", "lead", "water"]


def test_same_seed_reproducible_across_hash_seeds():
    """文字列ハッシュシードの異なる2プロセスで、同一seedの輸送結果が一致する。"""
    outs = []
    for hash_seed in ("0", "12345"):
        env = {**os.environ, "PYTHONHASHSEED": hash_seed}
        p = subprocess.run([sys.executable, "-c", _SCRIPT], capture_output=True,
                            text=True, env=env, check=True)
        outs.append(p.stdout)
    assert outs[0] == outs[1]
