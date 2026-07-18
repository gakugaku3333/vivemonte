# 光子エネルギー依存性EGS5相互検証結果（水、150 keV・高エネルギー側）

事前登録: [PREREGISTRATION.md](PREREGISTRATION.md)（結果を見る前に基準を固定済み）
関連: [docs/egs5_crosscheck/RESULTS.md](../RESULTS.md)（Phase 1, 60 keV水10cm）、
[water20kev/RESULTS.md](../water20kev/RESULTS.md)（対になる低エネルギー側）

## 実装前に確認した既知の落とし穴（xraylibスプライン外挿エラー）

`xraylib.CS_Photo_Partial`は軽元素・高エネルギー(E>約103keV)でスプライン
外挿エラーを起こすことがある。水は酸素(Z=8)を含むため、`chatcarlo/physics.py`
の`sample_fluorescence`が呼び出すこの関数が150 keVでエラーを起こさないか
事前に確認した。結果、酸素・水素とも K線最大エネルギーが蛍光カットオフ
(5 keV)未満のため既存のガード（`fluorescence_k_data`のmax<cutoffチェック）
で`CS_Photo_Partial`自体が呼ばれないことを確認し、ブロッカーにならないと
判断した（詳細はPREREGISTRATION.md参照）。実行結果もこの判断を裏付けている
（エラーなく完走）。

## 世界境界の空気路バイアスへの対応

[water20kev/RESULTS.md](../water20kev/RESULTS.md)で発覚した既知の落とし穴
（既定bbox_margin_cm=50cmによる隠れた空気減弱バイアス）を踏まえ、当初から
`Geometry(geoms, bbox_margin_cm=0.01)`とスラブ直近の線源配置で実行した。

## 物理モデル設定（両コード共通）

| 項目 | 設定 |
|---|---|
| 入射光子 | 単色150 keV、鉛筆ビーム |
| スラブ | 水、厚さ10 cm、密度1.0 g/cm³、真空境界 |
| history数 | 500,000（両コード） |
| seed / inseed | 1 |
| コンプトン散乱 | 束縛コンプトン（IBOUND=1） |
| レイリー散乱 | 有効（IRAYL=1） |
| PEGS5カットオフ | UE=0.700, UP=0.200 MeV（150 keVの光子・散乱電子を上回るよう引き上げ） |

## 結果: 一次透過率

| 手法 | 一次透過率 | 統計誤差 |
|---|---|---|
| 解析解（Beer-Lambert, xraylib由来μ=0.150524/cm） | 22.1964% | — |
| ChatCarlo MC（n=500,000, seed=1） | 22.2116% | ±0.0588pp（二項近似） |
| EGS5（RHO=1.0, n=500,000, inseed=1） | 22.26% | ±0.0588pp（二項近似） |

- ChatCarlo vs EGS5: 差0.0484pp、相対差**-0.217%**、**0.58σ**
- ChatCarlo vs 解析解: 相対差+0.068%
- EGS5 vs 解析解: 相対差+0.288%

## 事前基準との照合

**合格基準（事前登録）: 相互差2%以内かつ2σ以内。**

ChatCarlo-EGS5間の相対差0.217%・0.58σは基準を大きく下回る。**合格。**
Phase 1（60 keV、相対差0.88%・1.7σ）より一致度が高い。150 keVは
コンプトン散乱がさらに支配的になる領域だが、束縛コンプトン断面積モデル
（IBOUND=1で揃えた条件）は高エネルギー側でも良好に一致した。

## 判定

**合格。** 150 keV（診断領域上限）でChatCarloの一次透過率はEGS5と
統計誤差内で一致する。20 keV版（[water20kev/RESULTS.md](../water20kev/RESULTS.md)）
と合わせ、診断領域の下限・上限いずれでも、Phase 1（60 keV）で確立した
一致度と同等以上の結果が得られた。

## 監査結果

`vive-audit`スキル経由で4件合同のステージB監査を実施。**総合判定: 合格。**
一次データからの独立再計算で数値の一致を確認。要注意所見1件（PREREGISTRATION.md・
RESULTS.mdがgit未追跡でタイムスタンプのみが事前登録の根拠）を除き重大な指摘なし。
詳細は[energy_dependence_summary.md](../energy_dependence_summary.md)参照。

## 生データ・再現用ファイル

- ChatCarlo: [run_chatcarlo_water150.py](run_chatcarlo_water150.py) /
  [run_chatcarlo_water150.log](run_chatcarlo_water150.log)
- EGS5: [water150kev.f](water150kev.f) / [water150kev.inp](water150kev.inp) /
  [egs5job.out](egs5job.out) / [pgs5job.pegs5lst](pgs5job.pegs5lst)
