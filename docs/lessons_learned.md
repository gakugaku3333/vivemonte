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

## 線量係数（H*(10)等）は記憶で書かず1次資料を取得する

ICRP Publication 74 / ICRU Report 57 の h*(10)/Φ 表は数値を暗記していても
誤る危険が高い防護量データで、公式PDFの数値表は検索だけでは本文抽出できない
ことが多い（WebFetchがPDFのバイナリ内容を正しくテキスト化できないケースが
あった）。OpenMCプロジェクト（MITライセンス）がICRP74値をそのままテキスト
ファイルとして同梱していたため、`gh api repos/openmc-dev/openmc/contents/...`
で直接取得できた（[scripts/fetch_h_star_10.py](../scripts/fetch_h_star_10.py)）。

**対応**: 新しい線量換算係数を追加するときは、(1)記憶から値を書かない、
(2)1次資料PDFが読めない場合はオープンソースプロジェクトが同梱する
転記データ（ライセンス・出典明記の上で）を探す、(3)取得元をコード内
コメントとfetchスクリプトに明記する、の順で対応する。

## 「重元素ほど前方散乱に偏る」という直感は誤り（レイリー散乱の原子形状因子）

原子形状因子 F(Z,q) は q=0 で F=Z、単調減少するが、その減衰の速さは
**軽元素の方が急峻**（xraylibで実測: q=0.5 Å⁻¹でのF/ZはC:0.28, Ca:0.41, Pb:0.59）。
そのため同じ光子エネルギーでは軽元素の方が前方散乱に強く偏り、
重元素は後方まで広がる — 「重い原子ほど前方に強く散乱するはず」という
直感は逆だった。テストを先に直感で書いて失敗し、xraylibの生データを
確認してから修正した（[tests/test_rayleigh.py](../tests/test_rayleigh.py)）。

**対応**: 物理的な直感でテストのアサーションを書く前に、実装が使う
1次データ（この場合はxraylib自体）で数値を直接確認してから書く。
直感と実装のどちらが誤っているか区別がつくまでテストの向きを決め打ちしない。
