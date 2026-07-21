"""検証プロトコル(a): レイリー角度サンプリングの分布同等性フル検証。

docs/plan_rayleigh_compton_importance_sampling.md の検証プロトコル(a)に対応する
フル版（Phase 0/1で使うスクリプト。恒久回帰テストの縮小版は
tests/test_rayleigh_distribution.py）。

代表点マトリクス Z ∈ {1(H), 8(O), 20(Ca), 82(Pb)} × E ∈ {20, 60, 80, 150} keV の
全16組について:
  (1) 新実装（or Phase 0時点では現行実装のみ）が理論pdf
      p(cosθ) ∝ (1+cos²θ)F(Z,q)² に一標本適合度検定（KS）で適合するか
  (2) モーメント（⟨cosθ⟩, ⟨cos²θ⟩）が理論値と標準誤差3倍以内で一致するか
  (3)（新実装が存在する場合のみ）新旧2標本KS検定で同一分布とみなせるか

有意水準はBonferroni補正: α = 0.05/16 （16組の多重比較）。
「適合」の判定は p >= α（帰無仮説「理論分布/同一分布」を棄却しない）。

**サンプル数は新旧非対称**（2026-07-21、実測に基づく設計変更）: 当初は新旧とも
n=1,000,000を想定していたが、旧実装（cosθ一様提案）はZ=1(H)・E=150keVで平均
約4,872試行/光子（docs/speedup_baseline/RCI_PHASE0_BASELINE.md参照）に達し、
実測でn=100,000・この1組だけで390秒かかった——n=1,000,000なら1組だけで
約65分、4エネルギー×Hだけで数時間規模になり非現実的。旧実装を大きいnで
検証する意味も薄い（旧実装の低受理率そのものが本計画の置き換え理由であり、
高精度に検証すべき対象ではない）。そのため:
- 新実装（`sample_rayleigh_cos_theta`、Phase 1後）: n=1,000,000（高速、
  理論pdfとの一標本比較に十分な検出力）
- 旧実装（`_sample_rayleigh_cos_theta_uniform`、Phase 0時点の"current"、
  Phase 1後は"old"）: n=20,000（低受理率でも数十秒規模で完走する上限。
  2標本KS検定はscipy.stats.ks_2sampが不等サンプルサイズに対応しているため
  問題なく機能する。検出力は新実装側の一標本検定より低いが、Phase 0の目的
  はハーネス自体の較正であり、Phase 1の主たる合否判定は新実装vs理論の
  一標本検定が担う）

実行: .venv/bin/python docs/rci_verification/verify_rayleigh_distribution.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
from scipy import stats
from scipy.integrate import cumulative_trapezoid

from chatcarlo.materials import rayleigh_form_factor_table
from chatcarlo import physics

_HC_KEV_ANGSTROM = 12.3984193
_N_SAMPLE_FAST = 1_000_000  # 新実装（高受理率）用
_N_SAMPLE_SLOW = 20_000     # 旧実装（低受理率、cosθ一様提案）用
_N_GRID = 200_000
_Z_LIST = (1, 8, 20, 82)
_E_LIST = (20.0, 60.0, 80.0, 150.0)
_ALPHA = 0.05 / (len(_Z_LIST) * len(_E_LIST))  # Bonferroni補正
_MOMENT_SIGMA = 3.0


def _theoretical_cdf(z: int, e_keV: float):
    """理論分布 p(cosθ) ∝ (1+cos²θ)F(Z,q)² の累積分布関数を密な格子上で構築する。"""
    q_grid, f_grid = rayleigh_form_factor_table(z)
    c = np.linspace(-1.0, 1.0, _N_GRID)
    theta = np.arccos(c)
    q = e_keV * np.sin(theta / 2.0) / _HC_KEV_ANGSTROM
    f = np.interp(q, q_grid, f_grid)
    pdf_unnorm = (1.0 + c ** 2) * f ** 2
    cdf = cumulative_trapezoid(pdf_unnorm, c, initial=0.0)
    norm = cdf[-1]
    cdf /= norm
    pdf = pdf_unnorm / norm
    return c, cdf, pdf


def _theoretical_moments(c: np.ndarray, pdf: np.ndarray):
    mean_c = np.trapezoid(c * pdf, c)
    mean_c2 = np.trapezoid(c ** 2 * pdf, c)
    return mean_c, mean_c2


def _sample_moments_with_sem(samples: np.ndarray):
    n = len(samples)
    mean_c = samples.mean()
    sem_c = samples.std(ddof=1) / np.sqrt(n)
    c2 = samples ** 2
    mean_c2 = c2.mean()
    sem_c2 = c2.std(ddof=1) / np.sqrt(n)
    return mean_c, sem_c, mean_c2, sem_c2


def _get_samplers():
    """検証対象のサンプラーを返す: {"new": callable, "old": callable or None}。

    Phase 0時点（本ファイル作成時）では新実装はまだ存在しないため、現行実装
    (`sample_rayleigh_cos_theta`)のみを"current"として理論pdfと比較する。
    Phase 1で`_sample_rayleigh_cos_theta_uniform`（旧実装保持）が追加されたら
    自動的に新旧2標本比較も走る。
    """
    samplers = {"current": physics.sample_rayleigh_cos_theta}
    old = getattr(physics, "_sample_rayleigh_cos_theta_uniform", None)
    if old is not None and old is not physics.sample_rayleigh_cos_theta:
        samplers = {"new": physics.sample_rayleigh_cos_theta, "old": old}
    return samplers


def run():
    samplers = _get_samplers()
    mode = "Phase 1後（新旧比較モード）" if "new" in samplers else "Phase 0時点（現行実装単体モード）"
    print(f"=== レイリー角度サンプリング分布検証: {mode} ===")
    print(f"n_sample(fast/new)={_N_SAMPLE_FAST}, n_sample(slow/old)={_N_SAMPLE_SLOW}, "
          f"Bonferroni補正 alpha={_ALPHA:.6f}\n")

    all_pass = True
    rows = []
    for z in _Z_LIST:
        for e in _E_LIST:
            rng_seed = hash((z, int(e))) & 0xFFFFFFFF
            c_grid, cdf_grid, pdf_grid = _theoretical_cdf(z, e)
            theory_mean_c, theory_mean_c2 = _theoretical_moments(c_grid, pdf_grid)

            def cdf_fn(x, c_grid=c_grid, cdf_grid=cdf_grid):
                return np.interp(x, c_grid, cdf_grid)

            sample_results = {}
            for label, fn in samplers.items():
                # "current"はPhase 0時点では旧実装そのもの（低受理率）なのでSLOW側。
                # "new"（Phase 1後の高速新実装）だけがFAST側を使う。
                n_sample = _N_SAMPLE_FAST if label == "new" else _N_SAMPLE_SLOW
                rng = np.random.default_rng(rng_seed + (0 if label in ("current", "new") else 1))
                z_arr = np.full(n_sample, z)
                e_arr = np.full(n_sample, e)
                samples = fn(z_arr, e_arr, rng)
                ks = stats.kstest(samples, cdf_fn)
                mean_c, sem_c, mean_c2, sem_c2 = _sample_moments_with_sem(samples)
                moment_ok = (abs(mean_c - theory_mean_c) < _MOMENT_SIGMA * sem_c and
                             abs(mean_c2 - theory_mean_c2) < _MOMENT_SIGMA * sem_c2)
                fit_ok = ks.pvalue >= _ALPHA
                sample_results[label] = dict(samples=samples, ks=ks, mean_c=mean_c, sem_c=sem_c,
                                              mean_c2=mean_c2, sem_c2=sem_c2,
                                              fit_ok=fit_ok, moment_ok=moment_ok)
                status = "OK" if (fit_ok and moment_ok) else "NG"
                if not (fit_ok and moment_ok):
                    all_pass = False
                print(f"Z={z:>3} E={e:>5.0f}keV [{label:>7}] vs理論: KS p={ks.pvalue:.4f} "
                      f"({'fit_ok' if fit_ok else 'FIT_NG'}) "
                      f"<cosθ>={mean_c:+.4f}(理論{theory_mean_c:+.4f}) "
                      f"<cos2θ>={mean_c2:.4f}(理論{theory_mean_c2:.4f}) [{status}]")

            if "new" in sample_results and "old" in sample_results:
                ks2 = stats.ks_2samp(sample_results["new"]["samples"], sample_results["old"]["samples"])
                two_sample_ok = ks2.pvalue >= _ALPHA
                if not two_sample_ok:
                    all_pass = False
                print(f"Z={z:>3} E={e:>5.0f}keV [new vs old] 2標本KS p={ks2.pvalue:.4f} "
                      f"[{'OK' if two_sample_ok else 'NG'}]")
                rows.append((z, e, two_sample_ok))
            print()

    print("=== 総合結果:", "全通過" if all_pass else "不合格あり", "===")
    return all_pass


if __name__ == "__main__":
    ok = run()
    sys.exit(0 if ok else 1)
