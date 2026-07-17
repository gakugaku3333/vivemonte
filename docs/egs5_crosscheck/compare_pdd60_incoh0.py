"""Compare the standard EGS5 PDD/OCR run (INCOH=0, the established-code baseline,
docs/egs5_crosscheck/pdd60_phantom/egs5job.out) against the current ChatCarlo result
(docs/egs5_crosscheck/chatcarlo_pdd60_results.json). Mirrors compare_pdd60_incoh1.py's
methodology (mass_g = bin_volume_cm3 * RHO_EGS5=1.001, same KEV_PER_G_TO_GY constant
as the ChatCarlo-side script) but points at the standard (non-diagnostic) EGS5 run.
"""
import re, json, math

EGS5_OUT = "/Users/oishifamily/Projects/viveMonte/docs/egs5_crosscheck/pdd60_phantom/egs5job.out"
CHATCARLO_JSON = "/Users/oishifamily/Projects/viveMonte/docs/egs5_crosscheck/chatcarlo_pdd60_results.json"
KEV_PER_G_TO_GY = 1.602176634e-13
RHO_EGS5 = 1.001  # g/cm3, PEGS5 template density

with open(CHATCARLO_JSON) as f:
    cc = json.load(f)

pat = re.compile(
    r"^\s*(\S+)\s+mean\(MeV\)=\s*([\-0-9.EeD+]+)\s+sem\(MeV\)=\s*([\-0-9.EeD+]+)\s+relerr\(%\)=\s*([\-0-9.]+)"
)

rows = []
with open(EGS5_OUT) as f:
    for line in f:
        m = pat.match(line)
        if m:
            name, mean_s, sem_s, relerr_s = m.groups()
            rows.append((name, float(mean_s), float(sem_s), float(relerr_s)))

assert len(rows) == 47, f"expected 47 bins, got {len(rows)}"

def volume_cm3(name):
    # PDD bins: 2x2x1 cm = 4 cm3; lateral bins: 2x1x1 cm = 2 cm3
    if name.startswith("pdd_"):
        return 4.0
    else:
        return 2.0

results = []
for name, mean_mev, sem_mev, relerr in rows:
    vol = volume_cm3(name)
    mass_g = vol * RHO_EGS5
    egs5_gy = (mean_mev * 1000.0 / mass_g) * KEV_PER_G_TO_GY
    egs5_sem_gy = (sem_mev * 1000.0 / mass_g) * KEV_PER_G_TO_GY
    egs5_relerr_pct = 100.0 * egs5_sem_gy / egs5_gy

    cc_entry = cc[name]
    cc_gy = cc_entry["mean_Gy_per_history"]
    cc_sem_gy = cc_entry["sem_Gy_per_history"]
    cc_relerr_pct = 100.0 * cc_entry["rel_err"]

    reldiff_pct = 100.0 * (egs5_gy - cc_gy) / cc_gy
    sigma = math.sqrt(egs5_sem_gy**2 + cc_sem_gy**2)
    sigma_n = abs(egs5_gy - cc_gy) / sigma if sigma > 0 else float("nan")

    results.append({
        "name": name,
        "egs5_gy": egs5_gy,
        "egs5_relerr_pct": egs5_relerr_pct,
        "cc_gy": cc_gy,
        "cc_relerr_pct": cc_relerr_pct,
        "reldiff_pct": reldiff_pct,
        "sigma_n": sigma_n,
    })

print(f"{'bin':22s} {'EGS5(Gy/hist)':>14s} {'EGS5err%':>9s} {'CC(Gy/hist)':>14s} {'CCerr%':>8s} {'reldiff%':>9s} {'sigma':>7s}")
for r in results:
    print(f"{r['name']:22s} {r['egs5_gy']:.4e} {r['egs5_relerr_pct']:9.3f} {r['cc_gy']:.4e} {r['cc_relerr_pct']:8.3f} {r['reldiff_pct']:9.3f} {r['sigma_n']:7.3f}")

avg47 = sum(r["reldiff_pct"] for r in results) / len(results)
print()
print(f"47-bin average relative diff (EGS5-ChatCarlo)/ChatCarlo = {avg47:.4f} %")

# restrict to bins where both stat errors < 1% (same criterion as v1)
sub = [r for r in results if r["egs5_relerr_pct"] < 1.0 and r["cc_relerr_pct"] < 1.0]
print(f"bins with both stat err <1%: {len(sub)} / 47")
avg_sub = sum(r["reldiff_pct"] for r in sub) / len(sub)
print(f"avg reldiff over those bins = {avg_sub:.4f} %")

n_pass = sum(1 for r in sub if r["sigma_n"] < 2.0 and abs(r["reldiff_pct"]) < 2.0)
n_fail = sum(1 for r in sub if not (r["sigma_n"] < 2.0 and abs(r["reldiff_pct"]) < 2.0))
print(f"pre-registered criterion (sigma<2 AND |reldiff|<2%) among stat-err<1% bins: pass={n_pass} fail={n_fail}")

# PDD-only vs lateral-only subgroup averages
pdd_avg = sum(r["reldiff_pct"] for r in results if r["name"].startswith("pdd_")) / 15
print(f"PDD-only (15 bins) avg reldiff = {pdd_avg:.4f} %")
lat_avg = sum(r["reldiff_pct"] for r in results if not r["name"].startswith("pdd_")) / 32
print(f"lateral-only (32 bins) avg reldiff = {lat_avg:.4f} %")

minr = min(r["reldiff_pct"] for r in results)
maxr = max(r["reldiff_pct"] for r in results)
print(f"range: {minr:.3f}% to {maxr:.3f}%")

import csv
with open("/Users/oishifamily/Projects/viveMonte/docs/egs5_crosscheck/pdd60_phantom/comparison_table_pchip_boundcompton.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["bin","egs5_gy","egs5_err_pct","cc_gy","cc_err_pct","reldiff_pct","sigma"])
    for r in results:
        w.writerow([r["name"], r["egs5_gy"], r["egs5_relerr_pct"], r["cc_gy"], r["cc_relerr_pct"], r["reldiff_pct"], r["sigma_n"]])
