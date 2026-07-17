"""EGS5内部でINCOH=0とINCOH=1の差だけを直接比較する(ChatCarloを介さない)。

Phase 2b残差調査のvive-auditor監査(条件付き合格)の指摘事項に対応: 「対標準
EGS5(INCOH=0)で残る約-1.0%はEGS5側モデル簡略化に由来しChatCarlo側の課題では
ない」という結論を、過去フェーズの分析を借用せず、この2つのEGS5生ログのみから
独立に検算する。

もし ChatCarlo(束縛コンプトン)がEGS5(INCOH=1)と一致(+0.02%)しつつEGS5(INCOH=0)
と-1.05%ずれるなら、EGS5(INCOH=1)はEGS5(INCOH=0)よりおよそ+1.0〜1.1%高い値を
出しているはずである。これが確認できれば、残差がEGS5自身のINCOH設定の違いで
完全に説明できることが、ChatCarlo側の実装を一切参照せずに示せる。

実行: python3 compare_egs5_incoh0_vs_incoh1.py
"""
import re

EGS5_INCOH0 = "/Users/oishifamily/Projects/viveMonte/docs/egs5_crosscheck/pdd60_phantom/egs5job.out"
EGS5_INCOH1 = "/Users/oishifamily/Projects/viveMonte/docs/egs5_crosscheck/pdd60_phantom_incoh1/egs5job.out"

pat = re.compile(
    r"^\s*(\S+)\s+mean\(MeV\)=\s*([\-0-9.EeD+]+)\s+sem\(MeV\)=\s*([\-0-9.EeD+]+)\s+relerr\(%\)=\s*([\-0-9.]+)"
)


def read_rows(path):
    rows = {}
    with open(path) as f:
        for line in f:
            m = pat.match(line)
            if m:
                name, mean_s, sem_s, relerr_s = m.groups()
                rows[name] = (float(mean_s), float(sem_s), float(relerr_s))
    return rows


r0 = read_rows(EGS5_INCOH0)
r1 = read_rows(EGS5_INCOH1)
assert len(r0) == 47 and len(r1) == 47, (len(r0), len(r1))
assert set(r0) == set(r1)

print(f"{'bin':22s} {'INCOH0(MeV)':>14s} {'err0%':>7s} {'INCOH1(MeV)':>14s} {'err1%':>7s} {'reldiff%':>9s} {'sigma':>7s}")
diffs = []
diffs_lowerr = []
import math
for name in r0:
    m0, s0, e0 = r0[name]
    m1, s1, e1 = r1[name]
    reldiff = 100.0 * (m1 - m0) / m0
    sigma = math.sqrt(s0 ** 2 + s1 ** 2)
    sigma_n = abs(m1 - m0) / sigma if sigma > 0 else float("nan")
    diffs.append(reldiff)
    if e0 < 1.0 and e1 < 1.0:
        diffs_lowerr.append(reldiff)
    print(f"{name:22s} {m0:.4e} {e0:7.3f} {m1:.4e} {e1:7.3f} {reldiff:9.3f} {sigma_n:7.3f}")

print()
print(f"47-bin average (INCOH1-INCOH0)/INCOH0 = {sum(diffs)/len(diffs):.4f} %")
print(f"bins with both stat err <1%: {len(diffs_lowerr)}/47, avg = {sum(diffs_lowerr)/len(diffs_lowerr):.4f} %")
pdd_diffs = [d for n, d in zip(r0, diffs) if n.startswith("pdd_")]
lat_diffs = [d for n, d in zip(r0, diffs) if not n.startswith("pdd_")]
print(f"PDD-only avg = {sum(pdd_diffs)/len(pdd_diffs):.4f} %  (15 bins)")
print(f"lateral-only avg = {sum(lat_diffs)/len(lat_diffs):.4f} %  (32 bins)")
print()
print("参考(前回セッションの独立指標、ChatCarlo経由):")
print("  ChatCarlo(束縛)vsEGS5(INCOH=1) = +0.018%")
print("  ChatCarlo(束縛)vsEGS5(INCOH=0) = -1.050%")
implied = ((1 + 0.00018) / (1 - 0.01050) - 1) * 100
print(f"  上記2つから逆算されるEGS5(INCOH=1)/EGS5(INCOH=0)-1 の期待値 ≈ {implied:.3f}%")
