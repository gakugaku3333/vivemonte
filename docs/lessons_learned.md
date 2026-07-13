# 教訓

## xraylib の `CS_Energy`（μen/ρ）を吸収線量計算に使ってはいけない

`xraylib.CS_Energy` / `CS_Energy_CP` が返す質量エネルギー吸収係数 μen/ρ は、
NIST XAAMDI（Hubbell & Seltzer）の公表値と最大約17%（Al, 60 keVで実測）乖離する。
一方 μ/ρ（全断面積、`CS_Total`/`CS_Total_CP`）はNISTと1%以内で一致する。

**対応**: μ/ρ は xraylib、μen/ρ は NIST XAAMDI を直接取得したCSV
（`scripts/fetch_nist_xaamdi.py` → `vivemonte/data/nist_xaamdi/`）から
log-log補間で得る2系統構成にした（[materials.py](../vivemonte/materials.py)）。
`tests/test_materials.py` でNIST公表値とのスポット照合を自動テスト化済み。

**なぜ気づけたか**: 実装直後にNIST公表値と機械的に突き合わせるテストを書いたため。
新しい断面積・線量係数データソースを追加するときは、実装を信用する前に
必ず1次資料の数値と照合するテストを先に書くこと。

## 生体組織の組成規格の違いによるμ/ρの数%ずれは実在する物理差

xraylib の骨データは ICRP 組成、NIST XAAMDI は ICRU-44 組成で、
同じ「骨」でも μ/ρ が最大約2%異なる（60 keVで実測）。
これはバグではなく組成定義の違いによる正当な差なので、
テストの許容誤差は物質ごとに変える（`tests/test_materials.py` の `_LOOSE` 集合）。
