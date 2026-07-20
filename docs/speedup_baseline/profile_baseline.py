"""高速化計画(docs/plan_transport_speedup.md)の各Phase後に実行するプロファイル計測スクリプト。

実行: .venv/bin/python docs/speedup_baseline/profile_baseline.py <出力ファイル接頭辞>
      (リポジトリルートから。例: phase1 → phase1_timing.txt, phase1_profile.txt)
接頭辞省略時は baseline_*.txt に出力する（Phase 0専用。他Phaseで指定なしに実行すると
Phase 0の記録を上書きしてしまうため注意）。
"""
import cProfile
import io
import pstats
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

_PREFIX = sys.argv[1] if len(sys.argv) > 1 else "baseline"
_OUT_DIR = Path(__file__).resolve().parent

from chatcarlo.scene import load_scene
from chatcarlo.transport import run_transport

scene = load_scene("examples/chest_room.yaml")

# ウォームアップ（SpekPyキャッシュ・xraylibロード等の固定費を timing から除くため）
run_transport(scene, n_histories=1_000, seed=1)

timing_lines = []
for n in (1_000, 10_000, 100_000):
    t0 = time.time()
    run_transport(scene, n_histories=n, seed=42)
    dt = time.time() - t0
    line = f"n_histories={n:>8}  wall={dt:.3f}s  ({dt / n * 1e6:.2f} us/history)"
    print(line)
    timing_lines.append(line)

timing_path = _OUT_DIR / f"{_PREFIX}_timing.txt"
with open(timing_path, "w") as f:
    f.write("\n".join(timing_lines) + "\n")

pr = cProfile.Profile()
pr.enable()
run_transport(scene, n_histories=100_000, seed=42)
pr.disable()

s = io.StringIO()
ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
ps.print_stats(40)
profile_path = _OUT_DIR / f"{_PREFIX}_profile.txt"
with open(profile_path, "w") as f:
    f.write(s.getvalue())

print(f"\n--- saved: {timing_path}, {profile_path} ---")
