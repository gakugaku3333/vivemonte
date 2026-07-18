# K殻蛍光X線エネルギー依存性EGS5相互検証結果（銅、20 keV・低エネルギー側）

事前登録: [PREREGISTRATION.md](PREREGISTRATION.md)（結果を見る前に基準を固定済み）
関連: [docs/egs5_crosscheck/fluorescence_copper/RESULTS.md](../fluorescence_copper/RESULTS.md)
（100 keV、合格）

## 世界境界の空気路バイアスへの対応

[water20kev/RESULTS.md](../water20kev/RESULTS.md)で発覚した既知の落とし穴
（既定bbox_margin_cm=50cmによる低エネルギーでの隠れた空気減弱バイアス）を
踏まえ、`Geometry(geoms, bbox_margin_cm=0.01)`とスラブ直近の線源配置で
実行した（100 keV版の`fluorescence_copper/`は既定値のままだが、100 keVでは
影響が約1%と小さく問題化しなかった）。

## 物理モデル設定（両コード共通）

| 項目 | 設定 |
|---|---|
| 入射光子 | 単色20 keV、鉛筆ビーム |
| スラブ | 銅、厚さ0.01 cm（約3.03 mfp）、真空境界 |
| 密度 | 8.96 g/cm³（PEGS5組込みデフォルト8.9333と不一致のため明示指定、
  `pgs5job.pegs5lst`で確認済み） |
| history数 | 1,000,000（両コード） |
| seed / inseed | 1 |
| コンプトン散乱 | 束縛コンプトン（IBOUND=1） |
| レイリー散乱 | 有効（IRAYL=1） |
| Doppler broadening | 無効（ICPROF=0） |

## 観測量A: 脱出光子の全エネルギー透過率

統計誤差は標本標準偏差から直接計算（ChatCarlo:
`x.std(ddof=1)/√N`、EGS5: 0.1keV刻みヒストグラムから`E[X²]-E[X]²`で
分散を再構成）。

| ステップ | ChatCarlo | EGS5 | 相対差※ | σ数 | 事前基準 | 判定 |
|---|---|---|---|---|---|---|
| Step1（蛍光OFF） | 0.051608 ± 0.000221 | 0.051433 ± 0.000220 | +0.340% | 0.561σ | 2%以内かつ2σ以内 | **合格** |
| Step2（蛍光ON） | 0.073868 ± 0.000235 | 0.073986 ± 0.000234 | -0.160% | 0.356σ | 5%以内かつ3σ以内 | **合格** |

※相対差はEGS5を分母とする（(ChatCarlo−EGS5)/EGS5）。

### Step3（差分の一致、最も重視する）

- Δ_ChatCarlo = 0.073868 − 0.051608 = **0.022260**
- Δ_EGS5 = 0.073986 − 0.051433 = **0.022553**
- Δの比（ChatCarlo/EGS5） = **0.9870**

事前登録の基準（比0.7〜1.3倍以内）を明確に満たす。**合格。**

100 keV版（[fluorescence_copper/RESULTS.md](../fluorescence_copper/RESULTS.md)）では
Δがゼロ近傍で符号が逆転しフォールバック基準を使ったが、20 keVではΔが
両コードとも約2.2〜2.3%と大きく明確な値になり、比較としてより informativeで
明確な合格となった。20 keVでは入射光子の光電吸収分率自体が100 keVより
大幅に高く（銅のK端8.98 keVに近い）、K殻蛍光の生成・脱出への寄与が
統計的にはっきり検出できることが確認された。

## 観測量B: 脱出光子スペクトルの蛍光ピーク帯(7.8–9.2 keV)の割合

| 条件 | ChatCarlo | EGS5 |
|---|---|---|
| OFF | 0.000% | 0.008% |
| ON | 51.54% | 51.78% |

事前登録の予測（「自己吸収が強く観測量Bが非情報的になる」）は**外れた**。
20 keVでは薄いスラブ（0.01cm、20 keVでの3mfp）全体にわたって光電吸収が
分布するため、100 keV版（厚さ0.75cm、深部で生成された蛍光光子はほぼ
自己吸収される）とは異なり、蛍光ピーク帯の光子が実際に高い割合で脱出する。
蛍光ONでOFFの2倍以上という判定基準を大きく上回って満たし（ほぼ0%→約52%）、
かつ**ChatCarlo・EGS5の絶対値も0.24pp差と極めて近い**。**合格（判定に使用）。**

## 判定

**Step1: 合格、Step2: 合格、Step3（最重視）: 合格、観測量B: 合格。
総合判定: 合格。**

100 keV版の条件付き合格（銅は実際には全項目合格だった）よりもさらに
明確で、条件付きの要素がない完全な合格となった。事前登録の予測が外れた
点（観測量Bが非情報的にならなかった）は、20 keVという低エネルギー・
薄いスラブという組み合わせが100 keVとは異なる物理的レジーム（自己吸収より
生成分布の広がりが支配的）にあることを示しており、それ自体がエネルギー
依存性検証の主目的（異なるレジームでも実装が妥当か）に対する肯定的な結果。

## 監査結果

`vive-audit`スキル経由で4件合同のステージB監査を実施。**総合判定: 合格。**
一次データからの独立再計算で数値の一致を確認。要注意所見1件（PREREGISTRATION.md・
RESULTS.mdがgit未追跡でタイムスタンプのみが事前登録の根拠）を除き重大な指摘なし。
詳細は[energy_dependence_summary.md](../energy_dependence_summary.md)参照。

## 生データ・再現用ファイル

- ChatCarlo: [run_chatcarlo_fluorescence.py](run_chatcarlo_fluorescence.py) /
  `chatcarlo_results.json`
- EGS5: `egs5_off/fluor_off_cu20.f` / `egs5_off/fluor_off_cu20.inp` /
  `egs5_off/egs5job.out` / `egs5_off/pgs5job.pegs5lst`、
  `egs5_on/fluor_on_cu20.f` / `egs5_on/fluor_on_cu20.inp` /
  `egs5_on/egs5job.out` / `egs5_on/pgs5job.pegs5lst`
