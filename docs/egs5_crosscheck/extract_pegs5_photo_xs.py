"""PEGS5内部データから水の光電吸収断面積(60 keV)を直接抽出し、xraylibと比較する。

背景: PDD_RESULTS.md「原因の切り分け」節の仮説（ChatCarloのkerma推定量=NIST XAAMDI
μen/ρ参照 vs EGS5=PEGS5内部断面積による実MCシミュレーション、という2つのデータソース
の食い違い、特に光電吸収断面積の差がkerma/線量比較で増幅されている）を検証するため、
Phase 2a（BSF60_RESULTS.md「PEGS5断面積の直接照合」節、AIRDET空気の全断面積を
xraylibと相対0.02%で一致確認）と同じ手法を、水の光電吸収成分に適用する。

## 抽出方法

PEGS5(pegs5.f)は光子相互作用データを、対数エネルギー空間の区分線形フィット(PWLF)
として medium 生成物(pgs5job.pegs5lst / pgs5job.pegs5dat)に書き出す。photon側の
フィット対象は subroutine gfuns (pegs5.f 1716-1732行) が返す4量:

    v(1) = GMFP = 1/(pairtu+comptm+photte)      -- 全断面積(pair+compton+photo)の逆数
    v(2) = GBR1 = pair/(pair+compton+photo)     -- pair生成の累積確率
    v(3) = GBR2 = (pair+compton)/(pair+compton+photo)  -- pair+comptonの累積確率
    v(4) = (pair+compton+photo)/(pair+compton+photo+coherent) -- 非Rayleighである確率

これらは `call PWLF1(NGL,NALG,AP,UP,RMT2,EPG,ZTHRG,ZEPG,NIPG,DLOG,DEXP,AXG,BXG,
1000,4,AFG,BFG,GFUNS)` (pegs5.f 2972-2973行) で AFG(i,ifun),BFG(i,ifun) (ifun=1..4)
としてフィットされ、.pegs5lst に
  `$Echo write:EBINDA,BXG,AXG`                              (AXG,BXGは区分インデックス変換係数)
  `$Echo write:((bfg(i,ifun),afg(i,ifun),ifun=1,3),i=1,nge)` (GMFP,GBR1,GBR2)
  `$Echo write:(bfg(i,4),afg(i,4),i=1,nge)`                  (非Rayleigh確率)
としてダンプされる。

区分線形再構成 (QFIT, pegs5.f 4420-4515行): x=ln(E)として
  bin = int(AXG*x + BXG)
  value(x) = AFG(bin,ifun)*x + BFG(bin,ifun)

**重要な落とし穴**: comptm/photte/cohetm (pegs5.f) は `pcon=1.D-24*(AN*RHO/WM)*RLC`
という係数を使っており、これは断面積(barn→cm^2変換)×数密度に加えて **放射長RLC[cm]を
掛けている**。つまりPWLFがフィットしているtsansc=1/GMFPは「cm^-1」ではなく
「(cm^-1)×RLC」= 無次元量(放射長を単位とした不透明度)である。実際のマクロ断面積
Σ[cm^-1]を得るには `Σ = tsansc/RLC` としてRLCで割る必要がある(本モジュールの
主要な発見の一つ)。RLCは同じ.pegs5lstの `$ECHO WRITE:RLC,AE,AP,UE,UP` 行から取得。

## 検証(このモジュールでの自己チェック)

AIRDET(bsf60_free)の全断面積を同手法で再抽出し、Phase 2aが報告した0.18628 cm2/gと
比較 → 0.186277 cm2/g(相対差0.0002%)で再現。抽出手法・RLC補正が正しいことを確認。

実行: PYTHONPATH=. .venv/bin/python docs/egs5_crosscheck/extract_pegs5_photo_xs.py
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import xraylib

HERE = Path(__file__).resolve().parent
EGS5 = HERE / "egs5"


def _read_numbers_from(lines: list[str], start_idx: int, count: int) -> tuple[list[float], int]:
    nums: list[float] = []
    i = start_idx
    while len(nums) < count:
        for tok in lines[i].split():
            try:
                nums.append(float(tok))
            except ValueError:
                pass
        i += 1
    return nums[:count], i


def _find_line(lines: list[str], marker: str, start: int = 0) -> int:
    for i in range(start, len(lines)):
        if marker in lines[i]:
            return i
    raise ValueError(f"marker not found: {marker!r}")


def parse_pegs5lst(path: Path) -> dict:
    lines = path.read_text().splitlines()

    idx_nge = _find_line(lines, "NSGE,NGE,NSEKE")
    nge = int(lines[idx_nge + 1].split()[1])

    idx_rlc = _find_line(lines, "RLC,AE,AP,UE,UP")
    rlc = float(lines[idx_rlc + 1].split()[0])

    rho = None
    for line in lines:
        m = re.search(r"density=\s*([\d.]+)\s*\(g/cm\*\*3\)", line)
        if m:
            rho = float(m.group(1))
            break
    if rho is None:
        raise ValueError("density not found in pegs5lst")

    idx_axg = _find_line(lines, "EBINDA,BXG,AXG")
    (ebinda, bxg, axg), _ = _read_numbers_from(lines, idx_axg + 1, 3)

    idx_block13 = _find_line(lines, "bfg(i,ifun),afg(i,ifun),ifun=1,3")
    nums13, next_idx = _read_numbers_from(lines, idx_block13 + 1, nge * 3 * 2)
    afg = {ifun: [0.0] * (nge + 1) for ifun in (1, 2, 3, 4)}
    bfg = {ifun: [0.0] * (nge + 1) for ifun in (1, 2, 3, 4)}
    p = 0
    for i in range(1, nge + 1):
        for ifun in (1, 2, 3):
            bfg[ifun][i] = nums13[p]; p += 1
            afg[ifun][i] = nums13[p]; p += 1

    idx_block4 = _find_line(lines, "bfg(i,4),afg(i,4),i=1,nge", start=next_idx)
    nums4, _ = _read_numbers_from(lines, idx_block4 + 1, nge * 2)
    p = 0
    for i in range(1, nge + 1):
        bfg[4][i] = nums4[p]; p += 1
        afg[4][i] = nums4[p]; p += 1

    return dict(nge=nge, rlc=rlc, rho=rho, axg=axg, bxg=bxg, afg=afg, bfg=bfg)


def cross_sections_cm2_per_g(path: Path, energy_MeV: float) -> dict[str, float]:
    """PEGS5のPWLFフィットから、指定エネルギーでの成分別質量減弱係数[cm2/g]を再構成する。"""
    data = parse_pegs5lst(path)
    lnE = math.log(energy_MeV)
    idx = min(max(int(data["axg"] * lnE + data["bxg"]), 1), data["nge"])

    gmfp = data["afg"][1][idx] * lnE + data["bfg"][1][idx]  # 放射長単位
    gbr1 = data["afg"][2][idx] * lnE + data["bfg"][2][idx]
    gbr2 = data["afg"][3][idx] * lnE + data["bfg"][3][idx]
    nonrayl = data["afg"][4][idx] * lnE + data["bfg"][4][idx]

    tsansc_raw = 1.0 / gmfp  # = Sigma_macro[cm^-1] * RLC (無次元)
    cohr_raw = tsansc_raw * (1.0 / nonrayl - 1.0)

    rlc, rho = data["rlc"], data["rho"]
    sigma_pair = gbr1 * tsansc_raw / rlc
    sigma_compton = (gbr2 - gbr1) * tsansc_raw / rlc
    sigma_photo = (1.0 - gbr2) * tsansc_raw / rlc
    sigma_rayleigh = cohr_raw / rlc

    return {
        "pair": sigma_pair / rho,
        "compton": sigma_compton / rho,
        "photo": sigma_photo / rho,
        "rayleigh": sigma_rayleigh / rho,
        "total": (sigma_pair + sigma_compton + sigma_photo + sigma_rayleigh) / rho,
        "bin_index": idx,
    }


def main() -> None:
    E = 0.060  # MeV

    print("### 検証: AIRDET全断面積の再現(Phase 2aで報告された0.18628 cm2/gと比較) ###")
    airdet = cross_sections_cm2_per_g(EGS5 / "run_bsf60_free" / "pgs5job.pegs5lst", E)
    print(f"  PEGS5抽出 total = {airdet['total']:.6f} cm2/g  (Phase 2a報告値 0.18628, "
          f"相対差 {(airdet['total']-0.18628)/0.18628*100:+.4f}%)")
    print()

    print("### 水(H2O, RHO=1.001, IBOUND=1, 60 keV)の成分別質量減弱係数 ###")
    for name, sub in [("water60_bound (Phase 1)", "run_water60_bound"),
                       ("pdd60_phantom (Phase 2b)", "run_pdd60_phantom")]:
        path = EGS5 / sub / "pgs5job.pegs5lst"
        r = cross_sections_cm2_per_g(path, E)
        print(f"  [{name}] bin={r['bin_index']}")
        print(f"    photo={r['photo']:.6f}  compton={r['compton']:.6f}  "
              f"rayleigh={r['rayleigh']:.6f}  total={r['total']:.6f}  [cm2/g]")

    print()
    print("### xraylib (EPDL, CS_*_CP 'Water, Liquid', 60 keV) ###")
    photo_xr = xraylib.CS_Photo_CP("Water, Liquid", 60.0)
    compt_xr = xraylib.CS_Compt_CP("Water, Liquid", 60.0)
    rayl_xr = xraylib.CS_Rayl_CP("Water, Liquid", 60.0)
    total_xr = xraylib.CS_Total_CP("Water, Liquid", 60.0)
    print(f"    photo={photo_xr:.6f}  compton={compt_xr:.6f}  "
          f"rayleigh={rayl_xr:.6f}  total={total_xr:.6f}  [cm2/g]")

    r = cross_sections_cm2_per_g(EGS5 / "run_water60_bound" / "pgs5job.pegs5lst", E)
    print()
    print("### PEGS5 vs xraylib 相対差 (PEGS5-xraylib)/xraylib ###")
    for key, xr in [("photo", photo_xr), ("compton", compt_xr),
                     ("rayleigh", rayl_xr), ("total", total_xr)]:
        rel = (r[key] - xr) / xr * 100
        print(f"    {key:10s}: {rel:+.3f}%")


if __name__ == "__main__":
    main()
