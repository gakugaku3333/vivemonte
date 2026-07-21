"""Phase 0ベースライン: レイリー/コンプトン棄却ループの反復回数実測。

docs/speedup_baseline/RCI_PHASE0_BASELINE.md の代表点マトリクス表の再現用。
`_sample_rayleigh_cos_theta_uniform`/`sample_compton_bound`本体を呼ばず、
同じ棄却ロジックをここに複製して反復回数だけを数える（本体呼び出しだと
cos_theta配列の中身しか返らず、途中の反復回数が取れないため）。

実行: .venv/bin/python docs/rci_verification/measure_iters.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np

from chatcarlo.materials import incoherent_sq_table, rayleigh_form_factor_table

_HC_KEV_ANGSTROM = 12.3984193
_Z_LIST = (1, 8, 20, 82)
_E_LIST = (20, 60, 80, 150)


def count_rayleigh_iters(z: int, e_keV: float, n: int = 20_000, seed: int = 0):
    """旧実装(cosθ一様提案+包絡線2Z²棄却)の平均試行回数・ラウンド数を数える。"""
    rng = np.random.default_rng(seed)
    pending = np.arange(n)
    total_trials = 0
    rounds = 0
    q_grid, f_grid = rayleigh_form_factor_table(int(z))
    while len(pending) > 0:
        rounds += 1
        total_trials += len(pending)
        c = rng.uniform(-1.0, 1.0, len(pending))
        theta = np.arccos(c)
        q = e_keV * np.sin(theta / 2.0) / _HC_KEV_ANGSTROM
        f = np.interp(q, q_grid, f_grid)
        g = (1.0 + c ** 2) * f ** 2
        envelope = 2.0 * z ** 2
        xi2 = rng.random(len(pending))
        accept = xi2 * envelope <= g
        pending = pending[~accept]
    return total_trials / n, rounds


def count_compton_iters(z: int, e_keV: float, n: int = 20_000, seed: int = 0):
    """KN Kahn提案+S(Z,q)/Z追加棄却の平均試行回数・ラウンド数を数える。"""
    rng = np.random.default_rng(seed)
    mec2 = 511.0
    alpha = e_keV / mec2
    eps_min = 1.0 / (1.0 + 2.0 * alpha)
    envelope = 1.0 / eps_min + eps_min
    q_grid, s_grid = incoherent_sq_table(int(z))
    pending = np.arange(n)
    total_trials = 0
    rounds = 0
    while len(pending) > 0:
        rounds += 1
        total_trials += len(pending)
        xi1 = rng.random(len(pending))
        xi2 = rng.random(len(pending))
        eps_p = eps_min + xi1 * (1.0 - eps_min)
        cos_p = 1.0 - (1.0 / eps_p - 1.0) / alpha
        sin2_p = 1.0 - cos_p ** 2
        g = 1.0 / eps_p + eps_p - sin2_p
        accept_kn = xi2 * envelope <= g
        theta_p = np.arccos(np.clip(cos_p, -1.0, 1.0))
        q_p = e_keV * np.sin(theta_p / 2.0) / _HC_KEV_ANGSTROM
        s_over_z = np.interp(q_p, q_grid, s_grid) / z
        xi3 = rng.random(len(pending))
        accept = accept_kn & (xi3 <= s_over_z)
        pending = pending[~accept]
    return total_trials / n, rounds


if __name__ == "__main__":
    print("=== Rayleigh (旧実装: cosθ一様提案+包絡線2Z²) ===")
    for z in _Z_LIST:
        for e in _E_LIST:
            avg, rounds = count_rayleigh_iters(z, e)
            print(f"Z={z:>3} E={e:>3}keV: avg_trials={avg:>8.2f}  rounds={rounds}")

    print()
    print("=== Compton bound (KN Kahn + S/Z) ===")
    for z in _Z_LIST:
        for e in _E_LIST:
            avg, rounds = count_compton_iters(z, e)
            print(f"Z={z:>3} E={e:>3}keV: avg_trials={avg:>8.2f}  rounds={rounds}")
