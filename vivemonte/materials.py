"""材料データ — 断面積・減弱係数の取得。

- μ/ρ・内訳断面積（光電/コンプトン/レイリー）: xraylib（EPDLベース、NIST XCOMと一致確認済み）
  → 輸送カーネルの自由行程・相互作用抽選に使う
- μen/ρ: 同梱の NIST XAAMDI テーブル（scripts/fetch_nist_xaamdi.py で取得）
  → カーマ・吸収線量タリーに使う
  ※ xraylib の CS_Energy は NIST 公表値と最大約17%乖離するため使わない（検証済み）
"""
from __future__ import annotations

import functools
from pathlib import Path

import numpy as np
import xraylib

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


def _cs(func_elem, func_comp, material: str, energies_keV) -> np.ndarray:
    name, _, is_elem = resolve(material)
    e = np.atleast_1d(np.asarray(energies_keV, dtype=float))
    f = func_elem if is_elem else func_comp
    if is_elem:
        z = xraylib.SymbolToAtomicNumber(name)
        out = np.array([f(z, ek) for ek in e])
    else:
        out = np.array([f(name, ek) for ek in e])
    return out


def mu_rho(material: str, energies_keV) -> np.ndarray:
    """全質量減弱係数 μ/ρ [cm²/g]"""
    return _cs(xraylib.CS_Total, xraylib.CS_Total_CP, material, energies_keV)


@functools.lru_cache(maxsize=None)
def _load_xaamdi(key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = _DATA_DIR / f"{key}.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"NIST XAAMDIテーブルがありません: {path}\n"
            "scripts/fetch_nist_xaamdi.py を実行して取得してください")
    data = np.loadtxt(path, delimiter=",", comments="#")
    return data[:, 0], data[:, 1], data[:, 2]  # E_keV, mu/rho, muen/rho


def mu_en_rho(material: str, energies_keV) -> np.ndarray:
    """質量エネルギー吸収係数 μen/ρ [cm²/g]

    一次ソース: NIST XAAMDI（Hubbell & Seltzer）同梱テーブル、log-log補間。
    """
    key = _XAAMDI_FILES.get(material.lower().strip())
    if key is None:
        raise ValueError(
            f"材料 '{material}' の μen/ρ テーブルが未同梱です。"
            f"対応材料: {sorted(_XAAMDI_FILES)}。"
            "必要なら scripts/fetch_nist_xaamdi.py の TARGETS に追加してください")
    e_tab, _, muen_tab = _load_xaamdi(key)
    e = np.atleast_1d(np.asarray(energies_keV, dtype=float))
    if e.min() < e_tab[0] or e.max() > e_tab[-1]:
        raise ValueError(
            f"エネルギー {e.min():.3g}〜{e.max():.3g} keV はテーブル範囲 "
            f"[{e_tab[0]:.3g}, {e_tab[-1]:.3g}] keV 外です")
    return np.exp(np.interp(np.log(e), np.log(e_tab), np.log(muen_tab)))


def mu_rho_parts(material: str, energies_keV) -> dict[str, np.ndarray]:
    """内訳: 光電・コンプトン（非干渉性）・レイリー（干渉性） [cm²/g]"""
    return {
        "photoelectric": _cs(xraylib.CS_Photo, xraylib.CS_Photo_CP, material, energies_keV),
        "compton": _cs(xraylib.CS_Compt, xraylib.CS_Compt_CP, material, energies_keV),
        "rayleigh": _cs(xraylib.CS_Rayl, xraylib.CS_Rayl_CP, material, energies_keV),
    }


@functools.lru_cache(maxsize=None)
def element_composition(material: str) -> tuple[tuple[int, float], ...]:
    """材料 -> ((原子番号Z, 質量分率), ...)。単元素材料は1要素のタプル。

    レイリー散乱の角度分布は元素ごとの原子形状因子で決まるため、化合物・
    混合物ではどの構成元素で相互作用が起きたかを抽選する必要がある
    （transport.pyのレイリー散乱サンプリングで使用）。
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
    cs = np.array([[xraylib.CS_Rayl(int(z), ek) for ek in e] for z in zs])
    weighted = fracs[:, None] * cs
    total = weighted.sum(axis=0, keepdims=True)
    total = np.where(total > 0, total, 1.0)
    return zs, weighted / total


@functools.lru_cache(maxsize=None)
def rayleigh_form_factor_table(z: int, q_max: float = 20.0, n: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    """レイリー散乱の原子形状因子 F(Z,q) を q∈[0, q_max] Å⁻¹ でテーブル化（xraylib, EPDLベース）。

    q_max=20 Å⁻¹ は診断領域（kvp<=200keV、後方散乱θ=π）でも十分な余裕を持つ
    （E=200keV, θ=180°でも q≈16.1 Å⁻¹）。角度サンプリング側でnp.interpして使う。
    """
    q_grid = np.linspace(0.0, q_max, n)
    f_grid = np.array([xraylib.FF_Rayl(z, q) for q in q_grid])
    return q_grid, f_grid


def density(material: str) -> float:
    return resolve(material)[1]


def linear_mu(material: str, energies_keV) -> np.ndarray:
    """線減弱係数 μ [1/cm]"""
    return mu_rho(material, energies_keV) * density(material)
