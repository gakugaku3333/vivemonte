# 会話ワークフローによるPDD/OCR実演（学会抄録用、2026-07-16）

学会抄録（[../abstract_krmt2026.md](../abstract_krmt2026.md)）の【結果】記述——
「vive-interviewの会話ヒアリングでシーンを確定し、vive-checkの関門つき
ワークフローで計算する」——を、Phase 2b（[PDD_RESULTS.md](PDD_RESULTS.md)）の
水ファントムPDD/OCRで実際になぞった記録。目的は「会話だけで組んだ標準CLI
ワークフロー（`scene.yaml` + `chatcarlo run --dose-grid`）が、専用スクリプトで
既に監査済みのEGS5相互検証結果と同じ結論（形状一致）に行き着けるか」の実演。

## シーン

[../../examples/water_phantom_pdd_ocr.yaml](../../examples/water_phantom_pdd_ocr.yaml)。
Phase 2bと同条件（60 keV単色、10×10cm²、水ファントム30×20×30cm）だが、
Phase 2bの専用スクリプトが独自実装していた「平行ビーム近似」を、今回はscene.yaml
から直接指定できるよう`chatcarlo`本体に新機能として実装した:

- `source.spectrum: [{energy_keV: 60.0, weight: 1.0}]` — 単色ビーム
  （`source.kvp`の代替。両者は併用不可）
- `source.field.shape: parallel` — 非発散の平行ビーム（`size_cm`のみ、
  `sid_cm`不要。`rect`/`cone`と並ぶ第3の照射野形状）

両機能とも`mas`/`heel_effect`/`ctdi_vol_mGy`との併用は検証時にエラーになる
（それらはkvpベースのSpekPy計算に固定されており、上書きした条件と食い違うため）。
実装は[chatcarlo/scene.py](../../chatcarlo/scene.py)・
[chatcarlo/source.py](../../chatcarlo/source.py)、テストは
[tests/test_scene.py](../../tests/test_scene.py)・
[tests/test_parallel_beam.py](../../tests/test_parallel_beam.py)。

## ワークフロー

1. vive-interview: 目的（EGS5相互検証と同条件の実演）→撮影条件→ジオメトリー→
   計算条件（解像度1cm・n=3,000,000）を会話で確定。
2. vive-check 関門1: `validate`→`preview`→**vive-auditステージA監査**→
   ジオメトリー確認（ユーザー承認）。
3. vive-check 関門2: `trace`で光子軌跡確認（ユーザー承認）。
4. vive-check 関門3: 本計算（`-n 3000000 --dose-grid --resolution 1`）→
   **vive-auditステージB監査**。

## ステージB監査で発覚したバグと修正

初回の本計算はvive-auditorのステージB監査で**差し戻し**判定を受けた。
`--dose-grid`のtrack-lengthタリー（`accumulate_track_length`）が、平行ビームの
全光子が材料境界ちょうどから出発するという条件で、表面ボクセル層を約−2.7%
系統的に過小評価するバグを検出（詳細は
[../lessons_learned.md](../lessons_learned.md)「track-lengthタリーの
『サブステップ中点』スコアリングは…」節）。サブステップの中点スコアリングを
層化乱数点スコアリングに変更して修正し、**再監査で条件付き合格**。

このバグ発見・修正自体が、「関門つきワークフロー＋独立監査」という開発プロセスの
価値の実演にもなっている（監査なしでは表面層−2.7%のバイアスに気づかないまま
「EGS5と一致した」と誤報告するところだった——実際、修正前の初回計算では
たまたま2層目の値が表面ビン相当の値と偶然近く、見かけ上の一致に騙されかけた）。

## 結果（修正後、条件付き合格）

コマンド:
```
.venv/bin/python -m chatcarlo run examples/water_phantom_pdd_ocr.yaml \
    -n 3000000 --seed 42 --dose-grid --resolution 1 \
    --dose-out water_phantom_pdd_ocr_dose.npz
```

標準出力: 吸収0.5318・脱出0.4682・平均相互作用3.8989回/光子・
最大吸収線量4.78036e-15 Gy/history（表面層y=[0,1]cm、物理的に正しい位置）。

47ビン相当（PDD 15ビン＋OCR 2深さ×16ビン）に集約し、Phase 2b専用スクリプト・
EGS5と比較（[water_phantom_pdd_ocr_vs_egs5.png](water_phantom_pdd_ocr_vs_egs5.png)、
[water_phantom_pdd_ocr_maps.png](water_phantom_pdd_ocr_maps.png)）:

| 比較 | 内容 |
|---|---|
| 会話WF(CLI) vs Phase 2b専用スクリプト（ChatCarlo内部の2独立実装） | PDD平均差 −0.16%、最大差1.93%（統計誤差内） |
| 会話WF(CLI) vs EGS5、絶対値 | 平均 −1.60%（Phase 2bの既知の約−1.7%系統差を同水準で再現、未解明のまま） |
| 会話WF(CLI) vs EGS5、**形状**（各コード自身の基準ボクセルで正規化） | PDD平均+0.10%pt・最大0.44%pt／OCR表面平均+0.07%pt・最大1.25%pt／OCR 10cm平均+0.23%pt・最大1.19%pt |

**結論**: 会話ワークフローだけで組んだ標準CLI計算が、既監査済みのEGS5相互検証と
同じ結論（PDD/OCRの**形状**は一致、**絶対値**には未解明の約1.7%系統差が残る）に
独立に行き着いた。抄録の【結論】「確立コードEGS5との比較でも、水ファントム
深部・側方線量分布(PDD/OCR)の形状が2σの精度で一致し、物理的妥当性を確認した」は、
この実演で裏付けられている。

## 未解決・今後

- EGS5との絶対値−1.7%系統差の原因は依然未確定（[PDD_RESULTS.md](PDD_RESULTS.md)
  「原因の切り分け」参照、次点で優先度の高い調査項目）。
- 解像度0.5cm以下でのグリッド最大値成長の残存現象（新バグとは別機構、
  未検証）は[../lessons_learned.md](../lessons_learned.md)・
  [../../CLAUDE.md](../../CLAUDE.md)「Known sharp edges」節に記録済み。
