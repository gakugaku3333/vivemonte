"""材料データ — 断面積・減弱係数の取得。

- μ/ρ・内訳断面積（光電/コンプトン/レイリー）: xraylib（EPDLベース、NIST XCOMと一致確認済み）
  → 輸送カーネルの自由行程・相互作用抽選に使う。元素ごとにエネルギー格子上へ
  事前テーブル化しPCHIP補間で引く（docs/plan_transport_speedup.md Phase 1、
  下記「元素別断面積テーブル」参照）。混合物・化合物は質量分率で混合則
  （xraylibの*_CP関数と数値的に完全一致することを確認済み、下記参照）。
- μen/ρ: 同梱の NIST XAAMDI テーブル（scripts/fetch_nist_xaamdi.py で取得）
  → カーマ・吸収線量タリーに使う
  ※ xraylib の CS_Energy は NIST 公表値と最大約17%乖離するため使わない（検証済み）
"""
from __future__ import annotations

import functools
from pathlib import Path

import numpy as np
import xraylib
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

_DATA_DIR = Path(__file__).resolve().parent / "data" / "nist_xaamdi"

# 材料名 → 同梱NISTテーブルのファイル名（μen/ρ用）
_XAAMDI_FILES = {
    "water": "water", "air": "air", "soft_tissue": "soft_tissue",
    "bone": "bone", "lung": "lung", "muscle": "muscle", "adipose": "adipose",
    "pmma": "pmma", "concrete": "concrete", "lead_glass": "lead_glass",
    "aluminum": "z13_Al", "lead": "z82_Pb", "copper": "z29_Cu",
    "iron": "z26_Fe", "tungsten": "z74_W", "calcium": "z20_Ca",
}

# scene.yaml で使う短い材料名 → xraylib NIST化合物名
MATERIAL_ALIASES = {
    "water": "Water, Liquid",
    "air": "Air, Dry (near sea level)",
    "soft_tissue": "Tissue, Soft (ICRP)",
    "bone": "Bone, Cortical (ICRP)",
    "lung": "Lung (ICRP)",
    "pmma": "Polymethyl Methacralate (Lucite, Perspex)",
    "concrete": "Concrete, Portland",
    "lead": "Pb",
    "aluminum": "Al",
    "copper": "Cu",
    "iron": "Fe",
    "lead_glass": "Glass, Lead",
}

_DENSITY_OVERRIDE = {"Pb": 11.35, "Al": 2.699, "Cu": 8.96, "Fe": 7.874}


@functools.lru_cache(maxsize=None)
def resolve(material: str) -> tuple[str, float, bool]:
    """材料名 → (xraylib名, 密度 g/cm³, 元素かどうか)。未知ならValueError。"""
    name = MATERIAL_ALIASES.get(material.lower().strip(), material)
    if name in _DENSITY_OVERRIDE:
        return name, _DENSITY_OVERRIDE[name], True
    try:
        z = xraylib.SymbolToAtomicNumber(name)
        return name, xraylib.ElementDensity(z), True
    except ValueError:
        pass
    nist_names = xraylib.GetCompoundDataNISTList()
    if name in nist_names:
        return name, xraylib.GetCompoundDataNISTByName(name)["density"], False
    # あいまい一致で候補を提示（AIの自己修正用）
    cand = [n for n in nist_names if material.lower() in n.lower()][:5]
    raise ValueError(
        f"材料 '{material}' が見つかりません。候補: {cand or sorted(MATERIAL_ALIASES)}")


# --- 元素別断面積テーブル（docs/plan_transport_speedup.md Phase 1） ---
#
# 輸送カーネルは相互作用サンプリングのたびに μ/ρ・光電/コンプトン/レイリーの
# 内訳をxraylib経由で光子ごとに引いており、これがプロファイルで実測した
# 支配的コスト（chest_room.yaml n=1e5で全体4.5秒中 断面積呼び出し系だけで
# 約1.7秒、xraylib Cコール自体は数百万回）だった。μen/ρで既に確立済みの
# 「テーブル化＋PCHIP補間」パターンを断面積側にも横展開する。
#
# 元素(Z)単位でキャッシュする理由: 化合物のxraylib *_CP 関数は質量分率による
# 混合則（Σ w_i・cs_i(E)）と数値的に完全一致することを確認済み（Water/
# Concrete/Glass,Leadで検証、鉛ガラスのPb K端88keV直上でも一致）。元素単位で
# キャッシュすれば同じ元素を共有する材料間（水・骨・軟部組織のOなど）で
# テーブル構築を使い回せる。
#
# K/L/M...吸収端の扱い: 光電断面積は各殻の吸収端で不連続にジャンプする
# （物理的に正しい）。密なlog等間隔基本格子だけでは補間がこの不連続を
# なめらかに均してしまい歪める。そのためxraylibの全殻定数(_SHELL_CONSTANTS)
# を総当たりしてこの元素の全吸収端（K・L1-3・M1-5・...）を実測し、各端の
# 直下・直上に極小間隔の点を追加で挟むことで、補間区間そのものを不連続の
# 幅未満に狭め、ジャンプをほぼ段差のまま再現する。

_SHELL_CONSTANTS = [getattr(xraylib, _n) for _n in dir(xraylib) if _n.endswith("_SHELL")]

_XS_GRID_MIN_KEV = 1.0    # NIST XAAMDIテーブルの下限と揃える（mu_en_rhoと同じ制約）
_XS_GRID_MAX_KEV = 150.0  # 診断X線領域の上限（CLAUDE.md記載の設計スコープ「10–150 keV」
# と揃える）。scene.pyのkvp上限・source.spectrum明示指定の上限も同じ150 keVで検証する
# （2026-07-20、輸送カーネル高速化の副作用でテーブル化した際に判明: 以前はxraylib直呼び
# だったため上限が事実上無制限で、scene.pyのkvp検証は20〜200を許していた——これ自体が
# 「10–150 keV」という設計スコープと矛盾していた。テーブル化を機にスコープを150 keVに
# 統一し、範囲外はvalidate時にfail-fastで弾く設計に変更した）。
_XS_GRID_BASE_N = 2000    # log等間隔基本格子点数（隣接点間隔比約0.27%）
_EDGE_EPS_REL = 1e-5        # 実際のジャンプ位置を挟む点の相対オフセット（±0.001%）
_EDGE_SCAN_N = 200          # ジャンプ位置の粗スキャン点数
_EDGE_SCAN_CAP_REL = 0.02   # 粗スキャン窓の上限（±2%、隣接吸収端との最小間隔2.265%[Zn L2/L3]より狭い）
_EDGE_SCAN_MARGIN = 0.8     # 隣接吸収端までの距離に対する安全マージン（80%地点まで）
_EDGE_JUMP_RATIO_MIN = 1.2  # これ未満の比しか見つからなければ「ジャンプなし」とみなす
_NEAR_THRESHOLD_N = 30      # 端直後の近接閾値領域を補強する追加点数
_NEAR_THRESHOLD_CAP_REL = 0.05  # 補強範囲の上限（次の吸収端 or ジャンプ位置×5%の狭い方）


@functools.lru_cache(maxsize=None)
def _element_edge_energies_kev(z: int) -> tuple[float, ...]:
    """元素zの全吸収端の公称エネルギー（_XS_GRID_MIN_KEV〜_XS_GRID_MAX_KEV内のみ、
    xraylib.EdgeEnergy()の値。実際にCS_Photoがジャンプする位置とはわずかにズレる
    ため、格子点の配置には`_locate_photo_jump_kev`で実測した位置を使う）。
    """
    edges = []
    for shell in _SHELL_CONSTANTS:
        try:
            e = xraylib.EdgeEnergy(z, shell)
        except ValueError:
            continue
        if e and _XS_GRID_MIN_KEV < e < _XS_GRID_MAX_KEV:
            edges.append(float(e))
    return tuple(sorted(set(edges)))


def _locate_photo_jump_kev(z: int, nominal_edge_kev: float, window_rel: float) -> float:
    """CS_Photo(z,·)が実際にジャンプする位置を、公称吸収端の近傍で粗スキャン＋
    二分探索により特定する。

    実測で判明した注意点: xraylib.EdgeEnergy()の値と、実際にCS_Photoがジャンプする
    位置は完全には一致しない。固定の小さい探索窓（±0.2%）では不十分な元素があった
    （Zn L3端で実測約0.22%のズレを検出——EdgeEnergyでは1.0197keVだが実際の
    ジャンプは1.02〜1.0222keV付近）。そのため window_rel（呼び出し側が隣接吸収端
    との間隔から安全に決める）まで粗くスキャンしてジャンプ位置の見当をつけてから、
    その区間で二分探索により精密化する（docs/plan_transport_speedup.md参照）。
    ジャンプが見つからない場合は公称値にフォールバックする。
    """
    scan = np.geomspace(nominal_edge_kev * (1 - window_rel),
                         nominal_edge_kev * (1 + window_rel), _EDGE_SCAN_N)
    vals = np.array([xraylib.CS_Photo(z, e) for e in scan])
    ratios = vals[1:] / np.where(vals[:-1] > 0, vals[:-1], 1e-300)
    i = int(np.argmax(ratios))
    if ratios[i] < _EDGE_JUMP_RATIO_MIN:
        return nominal_edge_kev  # ジャンプが見つからない: 公称値にフォールバック
    a, b = scan[i], scan[i + 1]
    val_lo, val_hi = vals[i], vals[i + 1]
    mid_val = (val_lo + val_hi) / 2.0
    for _ in range(50):
        m = (a + b) / 2.0
        v = xraylib.CS_Photo(z, m)
        if v < mid_val:
            a = m
        else:
            b = m
    return (a + b) / 2.0


@functools.lru_cache(maxsize=None)
def _element_energy_grid_kev(z: int) -> np.ndarray:
    """元素zの断面積テーブル用エネルギー格子。基本log等間隔格子＋各吸収端の
    ジャンプ位置を挟む極小間隔の点に加え、ジャンプ直後の「近接閾値領域」を
    log等間隔で補強する。

    近接閾値領域の補強が必要な理由（実測で判明）: 亜鉛(Z=30)のL3端(1.02keV)と
    L2端(1.04keV)の間で、光電断面積はジャンプ後に単純な指数的減衰へは向かわず
    約2.3%のエネルギー幅でむしろ約1.8倍まで滑らかに上昇してから次の端へ至る
    （近接閾値構造。物理的に実在する挙動で、基本格子の密度でもこの区間内では
    捉えきれずPCHIP補間が最大約2%の誤差を生んだ）。ジャンプ直後〜次の吸収端
    （または相対5%上限）を追加点で密にすることで解消する。
    """
    nominal_edges = _element_edge_energies_kev(z)
    base = np.geomspace(_XS_GRID_MIN_KEV, _XS_GRID_MAX_KEV, _XS_GRID_BASE_N)
    pts = [base]
    for i, edge in enumerate(nominal_edges):
        # 探索窓は「固定上限2%」と「隣接吸収端までの距離の80%」の小さい方に
        # 制限し、別の吸収端を誤って跨がないようにする。
        neighbor_gaps = []
        if i > 0:
            neighbor_gaps.append((edge - nominal_edges[i - 1]) / edge)
        if i < len(nominal_edges) - 1:
            neighbor_gaps.append((nominal_edges[i + 1] - edge) / edge)
        window_rel = _EDGE_SCAN_CAP_REL
        if neighbor_gaps:
            window_rel = min(window_rel, min(neighbor_gaps) * _EDGE_SCAN_MARGIN)
        jump = _locate_photo_jump_kev(z, edge, window_rel)
        lo = max(_XS_GRID_MIN_KEV, jump * (1.0 - _EDGE_EPS_REL))
        hi = min(_XS_GRID_MAX_KEV, jump * (1.0 + _EDGE_EPS_REL))
        pts.append(np.array([lo, hi]))

        # 近接閾値領域の補強（ジャンプ直後 〜 次の吸収端 or +5%の狭い方）
        refine_cap = min(_XS_GRID_MAX_KEV, jump * (1.0 + _NEAR_THRESHOLD_CAP_REL))
        if i < len(nominal_edges) - 1:
            refine_cap = min(refine_cap, nominal_edges[i + 1])
        if refine_cap > hi:
            pts.append(np.geomspace(hi, refine_cap, _NEAR_THRESHOLD_N))
    return np.unique(np.concatenate(pts))


@functools.lru_cache(maxsize=None)
def _element_xs_tables(z: int) -> dict[str, np.ndarray]:
    """元素zの光電/コンプトン/レイリー質量減弱係数[cm²/g]を、格子上に
    一度だけxraylibで評価してlog-log空間の生テーブルとして返す（run内で使い回す）。

    lru_cacheでプロセス内キャッシュされるうえ、格子自体が密（隣接点間隔比
    約0.27%、吸収端近傍はさらに密）なので、区分線形補間で十分な精度が出る
    （tests/test_xs_tables.py参照。当初PchipInterpolatorで実装したが評価コストが
    xraylib呼び出し削減分を大きく相殺していた——chest_room.yaml n=1e5でPCHIP
    評価だけで約1.36秒——ため、素朴な補間に切り替えた）。3種の断面積が同じ
    エネルギー格子(log_e)を共有するので、後段で検索(searchsorted)を1回に
    まとめられるよう別々に格納する（docs/plan_transport_speedup.md Phase 1）。
    """
    grid = _element_energy_grid_kev(z)
    log_e = np.log(grid)
    photo = np.array([xraylib.CS_Photo(z, ek) for ek in grid])
    compt = np.array([xraylib.CS_Compt(z, ek) for ek in grid])
    rayl = np.array([xraylib.CS_Rayl(z, ek) for ek in grid])
    return {
        "log_e": log_e,
        "photo": np.log(photo),
        "compt": np.log(compt),
        "rayl": np.log(rayl),
    }


def _element_interp_index_frac(z: int, e: np.ndarray):
    """元素zの格子上での区分線形補間インデックス・重みを求める（光電/コンプトン/
    レイリーで共通のエネルギー格子を使うため、この検索(searchsorted)を1回だけ
    行い3種で使い回せば、np.interpを3回呼ぶより高速——同じ査問配列に対する
    探索を重複させないため）。
    """
    if e.min() < _XS_GRID_MIN_KEV or e.max() > _XS_GRID_MAX_KEV:
        raise ValueError(
            f"エネルギー {e.min():.3g}〜{e.max():.3g} keV は断面積テーブル範囲 "
            f"[{_XS_GRID_MIN_KEV:.3g}, {_XS_GRID_MAX_KEV:.3g}] keV 外です")
    log_e_grid = _element_xs_tables(z)["log_e"]
    log_e_query = np.log(e)
    idx = np.searchsorted(log_e_grid, log_e_query)
    idx = np.clip(idx, 1, len(log_e_grid) - 1)
    x0 = log_e_grid[idx - 1]
    x1 = log_e_grid[idx]
    frac = (log_e_query - x0) / (x1 - x0)
    return idx, frac


def _element_cs(z: int, energies_keV, kind: str) -> np.ndarray:
    """元素zの断面積[cm²/g]をテーブル補間で引く。kind: 'photo'/'compt'/'rayl'。"""
    e = np.atleast_1d(np.asarray(energies_keV, dtype=float))
    idx, frac = _element_interp_index_frac(z, e)
    y = _element_xs_tables(z)[kind]
    return np.exp(y[idx - 1] + frac * (y[idx] - y[idx - 1]))


def _element_cs_all(z: int, energies_keV) -> dict[str, np.ndarray]:
    """元素zの光電/コンプトン/レイリー断面積[cm²/g]を1回の探索でまとめて求める。"""
    e = np.atleast_1d(np.asarray(energies_keV, dtype=float))
    idx, frac = _element_interp_index_frac(z, e)
    tables = _element_xs_tables(z)
    return {kind: np.exp(tables[kind][idx - 1] + frac * (tables[kind][idx] - tables[kind][idx - 1]))
            for kind in ("photo", "compt", "rayl")}


def _compound_cs_all(material: str, energies_keV) -> dict[str, np.ndarray]:
    """化合物・混合物の光電/コンプトン/レイリー断面積[cm²/g]を、構成元素の
    質量分率加重和（混合則）でまとめて求める。"""
    comp = element_composition(material)
    e = np.atleast_1d(np.asarray(energies_keV, dtype=float))
    total = {"photo": np.zeros_like(e), "compt": np.zeros_like(e), "rayl": np.zeros_like(e)}
    for z, frac_mass in comp:
        parts = _element_cs_all(z, e)
        for kind in total:
            total[kind] = total[kind] + frac_mass * parts[kind]
    return total


def _cs_parts(material: str, energies_keV) -> dict[str, np.ndarray]:
    name, _, is_elem = resolve(material)
    if is_elem:
        z = xraylib.SymbolToAtomicNumber(name)
        parts = _element_cs_all(z, energies_keV)
    else:
        parts = _compound_cs_all(material, energies_keV)
    return {
        "photoelectric": parts["photo"],
        "compton": parts["compt"],
        "rayleigh": parts["rayl"],
    }


def mu_rho(material: str, energies_keV) -> np.ndarray:
    """全質量減弱係数 μ/ρ [cm²/g]"""
    parts = _cs_parts(material, energies_keV)
    return parts["photoelectric"] + parts["compton"] + parts["rayleigh"]


@functools.lru_cache(maxsize=None)
def _load_xaamdi(key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = _DATA_DIR / f"{key}.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"NIST XAAMDIテーブルがありません: {path}\n"
            "scripts/fetch_nist_xaamdi.py を実行して取得してください")
    data = np.loadtxt(path, delimiter=",", comments="#")
    return data[:, 0], data[:, 1], data[:, 2]  # E_keV, mu/rho, muen/rho


@functools.lru_cache(maxsize=None)
def _muen_interpolator(key: str) -> PchipInterpolator:
    e_tab, _, muen_tab = _load_xaamdi(key)
    return PchipInterpolator(np.log(e_tab), np.log(muen_tab))


def mu_en_rho(material: str, energies_keV) -> np.ndarray:
    """質量エネルギー吸収係数 μen/ρ [cm²/g]

    一次ソース: NIST XAAMDI（Hubbell & Seltzer）同梱テーブル、log-log PCHIP補間
    （格子点を厳密に通過する形状保存補間。log-log線形補間は30-80keV帯の格子
    間隔20keV区間で最大約3.3%の曲率誤差を持つことが判明したため採用、
    docs/egs5_crosscheck/pdd60_NOTES.md参照）。
    """
    key = _XAAMDI_FILES.get(material.lower().strip())
    if key is None:
        raise ValueError(
            f"材料 '{material}' の μen/ρ テーブルが未同梱です。"
            f"対応材料: {sorted(_XAAMDI_FILES)}。"
            "必要なら scripts/fetch_nist_xaamdi.py の TARGETS に追加してください")
    e_tab, _, _ = _load_xaamdi(key)
    e = np.atleast_1d(np.asarray(energies_keV, dtype=float))
    if e.min() < e_tab[0] or e.max() > e_tab[-1]:
        raise ValueError(
            f"エネルギー {e.min():.3g}〜{e.max():.3g} keV はテーブル範囲 "
            f"[{e_tab[0]:.3g}, {e_tab[-1]:.3g}] keV 外です")
    return np.exp(_muen_interpolator(key)(np.log(e)))


def mu_rho_parts(material: str, energies_keV) -> dict[str, np.ndarray]:
    """内訳: 光電・コンプトン（非干渉性）・レイリー（干渉性） [cm²/g]"""
    return _cs_parts(material, energies_keV)


@functools.lru_cache(maxsize=None)
def element_composition(material: str) -> tuple[tuple[int, float], ...]:
    """材料 -> ((原子番号Z, 質量分率), ...)。単元素材料は1要素のタプル。

    レイリー散乱の角度分布は元素ごとの原子形状因子で決まるため、化合物・
    混合物ではどの構成元素で相互作用が起きたかを抽選する必要がある
    （physics.pyのレイリー散乱サンプリングで使用）。
    """
    name, _, is_elem = resolve(material)
    if is_elem:
        return ((xraylib.SymbolToAtomicNumber(name), 1.0),)
    data = xraylib.GetCompoundDataNISTByName(name)
    return tuple(zip(data["Elements"], data["massFractions"]))


def rayleigh_element_weights(material: str, energies_keV) -> tuple[np.ndarray, np.ndarray]:
    """材料内でレイリー相互作用がどの構成元素で起きたかの重み。

    戻り値: (Z配列(n_elem,), 重み行列(n_elem, n_energies))。各列(energies_keVの
    1点ごと)の和が1になるよう、質量分率×元素別レイリー断面積で規格化する。
    """
    comp = element_composition(material)
    zs = np.array([z for z, _ in comp])
    fracs = np.array([f for _, f in comp])
    e = np.atleast_1d(np.asarray(energies_keV, dtype=float))
    cs = np.array([_element_cs(int(z), e, "rayl") for z in zs])
    weighted = fracs[:, None] * cs
    total = weighted.sum(axis=0, keepdims=True)
    total = np.where(total > 0, total, 1.0)
    return zs, weighted / total


@functools.lru_cache(maxsize=None)
def rayleigh_form_factor_table(z: int, q_max: float = 20.0, n: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    """レイリー散乱の原子形状因子 F(Z,q) を q∈[0, q_max] Å⁻¹ でテーブル化（xraylib, EPDLベース）。

    q_max=20 Å⁻¹ は診断領域（kvp<=150keV、後方散乱θ=π）でも十分な余裕を持つ
    （E=150keV, θ=180°でも q≈12.1 Å⁻¹）。角度サンプリング側でnp.interpして使う。
    """
    q_grid = np.linspace(0.0, q_max, n)
    f_grid = np.array([xraylib.FF_Rayl(z, q) for q in q_grid])
    return q_grid, f_grid


@functools.lru_cache(maxsize=None)
def rayleigh_cumulative_table(z: int, q_max: float = 20.0, n: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    """レイリー角度サンプリングの逆変換用累積分布テーブル（docs/
    plan_rayleigh_compton_importance_sampling.md の2段階方式、Phase 1）。

    変数変換 x ≡ q² を使うと cosθ = 1 - 2x/x_max(E) がxの1次式になり
    （ヤコビアン|dcosθ/dx|が定数）、目標分布 (1+cos²θ)F(Z,q)² は
    「F(Z,√x)²からの逆変換」と「角度因子(1+cos²θ)/2の棄却」に分解できる
    （`physics.sample_rayleigh_cos_theta`参照）。ここではその第1段用に
    A(x) = ∫₀ˣ F(Z,√x')² dx' を x∈[0, q_max²]上に事前計算する。

    `rayleigh_form_factor_table`と同じq_grid(線形, [0,q_max])をx=q²に写して
    使う——q間隔が一様でもx間隔はΔx≈2qΔqとqに比例して広がるため、x空間では
    小x(前方散乱・q小)側が自動的に密になり、F²の変化が速い領域をちょうど
    捉えられる。
    """
    q_grid, f_grid = rayleigh_form_factor_table(z, q_max=q_max, n=n)
    x_grid = q_grid ** 2
    a_grid = cumulative_trapezoid(f_grid ** 2, x_grid, initial=0.0)
    return x_grid, a_grid


def compton_element_weights(material: str, energies_keV) -> tuple[np.ndarray, np.ndarray]:
    """材料内でコンプトン相互作用がどの構成元素で起きたかの重み。

    `rayleigh_element_weights`と同型。質量分率×元素別コンプトン断面積
    （xraylib.CS_Compt、束縛効果込み）で規格化する。
    """
    comp = element_composition(material)
    zs = np.array([z for z, _ in comp])
    fracs = np.array([f for _, f in comp])
    e = np.atleast_1d(np.asarray(energies_keV, dtype=float))
    cs = np.array([_element_cs(int(z), e, "compt") for z in zs])
    weighted = fracs[:, None] * cs
    total = weighted.sum(axis=0, keepdims=True)
    total = np.where(total > 0, total, 1.0)
    return zs, weighted / total


def photo_element_weights(material: str, energies_keV) -> tuple[np.ndarray, np.ndarray]:
    """材料内で光電相互作用がどの構成元素で起きたかの重み（蛍光X線サンプリング用）。

    `compton_element_weights`/`rayleigh_element_weights`と同型。質量分率×
    元素別光電断面積（xraylib.CS_Photo）で規格化する。
    """
    comp = element_composition(material)
    zs = np.array([z for z, _ in comp])
    fracs = np.array([f for _, f in comp])
    e = np.atleast_1d(np.asarray(energies_keV, dtype=float))
    cs = np.array([_element_cs(int(z), e, "photo") for z in zs])
    weighted = fracs[:, None] * cs
    total = weighted.sum(axis=0, keepdims=True)
    total = np.where(total > 0, total, 1.0)
    return zs, weighted / total


@functools.lru_cache(maxsize=None)
def fluorescence_k_data(z: int) -> tuple[float, float, np.ndarray, np.ndarray]:
    """K殻蛍光データ（xraylib、EPDLベース）。

    戻り値: (K吸収端 keV, K蛍光収率 ω_K, 線エネルギー配列 keV, 線発生確率配列)
    線発生確率は8線（Kα2/Kα1: KL2/KL3、Kβ3/Kβ1: KM2/KM3、Kβ2系列: KN2/KN3、
    Kβ4・Kβ5相当の高殻由来分をまとめたKO_LINE/KP_LINE集計値）のRadRateで
    規格化した和=1の確率。全元素でK放射遷移の99.99%以上をカバーする（xraylibで
    Pb/W/Fe/Cu/Ca/Uにつき実測確認済み。当初KL2/KL3/KM2/KM3の4線のみだったが、
    Pbで95.05%しかカバーしておらず、EGS5相互検証の脱出光子スペクトルに
    88 keV付近の未説明ピークとして現れたことで発覚した欠落
    — docs/egs5_crosscheck/fluorescence/RESULTS.md参照）。
    線が存在しない（エラーまたは0を返す）場合はスキップし、有効線が1本も
    なければ ω_K=0.0 として扱う（軽元素・K蛍光を実質無視できる場合）。
    """
    edge_keV = xraylib.EdgeEnergy(z, xraylib.K_SHELL)
    lines = [xraylib.KL2_LINE, xraylib.KL3_LINE, xraylib.KM2_LINE, xraylib.KM3_LINE,
             xraylib.KN2_LINE, xraylib.KN3_LINE, xraylib.KO_LINE, xraylib.KP_LINE]
    energies = []
    rates = []
    for line in lines:
        try:
            rate = xraylib.RadRate(z, line)
            e_line = xraylib.LineEnergy(z, line)
        except ValueError:
            continue
        if rate > 0 and e_line > 0:
            energies.append(e_line)
            rates.append(rate)
    if not rates:
        return edge_keV, 0.0, np.array([]), np.array([])
    rates_arr = np.array(rates)
    rates_arr = rates_arr / rates_arr.sum()
    try:
        omega_k = xraylib.FluorYield(z, xraylib.K_SHELL)
    except ValueError:
        omega_k = 0.0
    return edge_keV, omega_k, np.array(energies), rates_arr


@functools.lru_cache(maxsize=None)
def incoherent_sq_table(z: int, q_max: float = 20.0, n: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    """コンプトン散乱の非干渉性散乱関数 S(Z,q) を q∈[0, q_max] Å⁻¹ でテーブル化
    （xraylib.SF_Compt、EPDLベース）。

    `rayleigh_form_factor_table`と対称的な用途: S(Z,q)はq→∞でZに単調収束し、
    束縛コンプトン散乱の追加棄却（S(Z,q)/Zを受理確率とする）に使う。xraylibの
    SF_Comptはq<~1e-3付近でスプライン外挿エラーになるため、q=0(物理的極限
    S=0)を手動で追加し、テーブルはq=1e-3から張る。
    """
    q_grid = np.concatenate([[0.0], np.linspace(1e-3, q_max, n - 1)])
    s_grid = np.concatenate([[0.0],
                              [xraylib.SF_Compt(z, q) for q in q_grid[1:]]])
    return q_grid, s_grid


def density(material: str) -> float:
    return resolve(material)[1]


def material_groups(names: np.ndarray):
    """材料名配列を (材料名, ブールマスク) の組に分けて材料名順に返す。

    輸送・タリーでは光子バッチを材料ごとにまとめて断面積を引く処理が
    頻出するため、そのグループ化を一箇所に集約する。ソートは必須:
    set()の反復順は文字列ハッシュのランダム化でプロセスごとに変わるため、
    グループ順に依存して乱数を消費する処理（レイリー元素抽選など）が
    未ソートだと同一seedでも実行ごとに結果が変わってしまう。
    """
    for name in sorted(set(names.tolist())):
        yield name, names == name


def linear_mu(material: str, energies_keV) -> np.ndarray:
    """線減弱係数 μ [1/cm]"""
    return mu_rho(material, energies_keV) * density(material)
