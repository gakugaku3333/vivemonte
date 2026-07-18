# 光子エネルギー依存性EGS5相互検証結果（水、20 keV・低エネルギー側）

事前登録: [PREREGISTRATION.md](PREREGISTRATION.md)（結果を見る前に基準を固定済み）
関連: [docs/egs5_crosscheck/RESULTS.md](../RESULTS.md)（Phase 1, 60 keV水10cm）、
[water150kev/RESULTS.md](../water150kev/RESULTS.md)（対になる高エネルギー側）

## 実装作業中に発見した落とし穴: 世界境界の空気路による隠れた減弱バイアス

ChatCarlo側スクリプトを既存60 keV版と同じ書き方（`Geometry`のbbox_margin_cm
既定50cm、線源をスラブ手前10cmに配置）で作ったところ、Beer-Lambert解析解
29.6787%に対しChatCarlo MC実測が28.08%と、統計誤差(±0.06pp)の25σ相当という
説明不能な乖離が出た。原因は、光子がスラブを抜けた後も既定bbox_marginにより
世界境界までさらに約50cmの空気中を飛行しており、20 keVではこの空気路減弱
（60cm合計で約5.5%）が無視できないため。`Geometry(geoms, bbox_margin_cm=0.01)`
とスラブ直近の線源配置に変更し、EGS5のtutor5パターン（スラブ外は真空）に
揃えたところ、解析解と統計誤差以内で一致した。**この教訓は
[docs/lessons_learned.md](../../lessons_learned.md)に別途記録し、60 keV版
Phase 1・100 keV蛍光検証群にも同種の（より小さい）バイアスが存在していた
可能性を指摘した。**

## 物理モデル設定（両コード共通）

| 項目 | 設定 |
|---|---|
| 入射光子 | 単色20 keV、鉛筆ビーム |
| スラブ | 水、厚さ1.5 cm、密度1.0 g/cm³、真空境界 |
| history数 | 500,000（両コード） |
| seed / inseed | 1 |
| コンプトン散乱 | 束縛コンプトン（IBOUND=1） |
| レイリー散乱 | 有効（IRAYL=1） |

## 結果: 一次透過率

| 手法 | 一次透過率 | 統計誤差 |
|---|---|---|
| 解析解（Beer-Lambert, xraylib由来μ=0.809828/cm） | 29.6787% | — |
| ChatCarlo MC（n=500,000, seed=1） | 29.7122% | ±0.0646pp（二項近似） |
| EGS5（RHO=1.0, n=500,000, inseed=1） | 29.76% | ±0.0647pp（二項近似） |

- ChatCarlo vs EGS5: 差0.0478pp、相対差**-0.161%**、**0.52σ**
- ChatCarlo vs 解析解: 相対差+0.113%
- EGS5 vs 解析解: 相対差+0.274%

## 事前基準との照合

**合格基準（事前登録）: 相互差2%以内かつ2σ以内。**

ChatCarlo-EGS5間の相対差0.161%・0.52σは基準を大きく下回る。**合格。**
Phase 1（60 keV、相対差0.88%・1.7σ）より一致度が高い。20 keVは光電効果が
より支配的になる領域（水の光電分率は60 keVより高い）だが、その領域でも
EGS5との一致は良好であり、事前予測（「IBOUND由来のコンプトン差の寄与が
相対的に小さくなる分、光電断面積自体のライブラリ差が相対的に見えやすくなる
可能性がある」）は今回のデータでは顕在化しなかった。

## 判定

**合格。** 20 keV（診断領域下限付近）でChatCarloの一次透過率はEGS5と
統計誤差内で一致する。

## 監査結果

`vive-audit`スキル経由で4件合同のステージB監査を実施。**総合判定: 合格。**
一次データからの独立再計算で数値の一致を確認。要注意所見1件（PREREGISTRATION.md・
RESULTS.mdがgit未追跡でタイムスタンプのみが事前登録の根拠）を除き重大な指摘なし。
詳細は[energy_dependence_summary.md](../energy_dependence_summary.md)参照。

## 生データ・再現用ファイル

- ChatCarlo: [run_chatcarlo_water20.py](run_chatcarlo_water20.py) /
  [run_chatcarlo_water20.log](run_chatcarlo_water20.log)
- EGS5: [water20kev.f](water20kev.f) / [water20kev.inp](water20kev.inp) /
  [egs5job.out](egs5job.out) / [pgs5job.pegs5lst](pgs5job.pegs5lst)
  （EGS5本体は第三者コードのためリポジトリに含めない、`.gitignore`参照）
