# 計画: 人間と確認しながら進める視覚ワークフロー（trace / plot / vive-checkスキル）

作成日: 2026-07-14 / ステータス: **計画のみ（未着手）** / 実装担当: Sonnet想定

## 目的

「AIが組んだ体系を人間が目視確認してから先へ進む」というviveMonteの
中核思想を、ジオメトリー確認だけでなく計算パイプライン全体に広げる。
3つの関門を作る:

1. **ジオメトリー確認** — 既存の `vivemonte preview`（実装済み、変更不要）
2. **光子軌跡確認** — 小historyで軌跡を3D表示する `vivemonte trace`（新規）
3. **結果確認** — 線量/H*(10)マップを断面図にする `vivemonte plot`（新規）

最後に、この3関門を順に通す進行役スキル `/vive-check` を作る。

## 設計原則（実装前に必ず読むこと）

- **能力はリポジトリに、進行だけスキルに**。trace/plotはCLIサブコマンド＋
  テスト付きのコードとして実装する。スキル本文には実装を持たせない
  （スキルはテストできず品質保証の網から漏れるため）。
- **既存の輸送カーネルの性能を劣化させない**。軌跡記録はopt-in
  （引数未指定時は完全に従来動作。既存55テストが無変更で通ること）。
- **preview.pyのHTMLテンプレートを複製しない**。軌跡表示は既存
  `_TEMPLATE` の拡張として実装し、previewとtraceで同一テンプレートを共有する
  （`DATA.trajectories` が空/未定義ならpreviewと同じ表示になる設計）。
- 新しい物理は一切追加しない。可視化と進行のみ。

---

## Deliverable 1: `vivemonte trace` — 光子軌跡の3D可視化

### CLI仕様

```bash
python -m vivemonte trace scene.yaml [-n 200] [--seed 42] [-o out.html]
```

- `-n` 既定200。**2000超は警告を出して2000にクランプ**する
  （HTML内蔵JSONの肥大とcanvas描画性能のため。1光子あたり平均3〜4区間
  ×2端点×3座標で、2000光子なら区間数万のオーダーに収まる）。
- `-o` 未指定時は `<scene>_trace.html`。
- 終了時に「光子N個の軌跡を書き出しました: out.html」と表示。

### 実装 1-a: 軌跡レコーダ（vivemonte/transport.py）

`transport_photons` に opt-in 引数を追加する:

```python
@dataclass
class TrajectoryRecorder:
    """軌跡記録（小history可視化用）。ループ1周ごとに飛行区間を追記する。

    segments: list of (start(3,), end(3,), energy_keV, event) を
    numpy配列のリストとして貯め、finalize()で結合する。
    event は区間の終端で起きたことを表す文字列:
      "boundary"（材料境界を通過して継続）, "photoelectric", "compton",
      "rayleigh", "escape"
    """
    starts: list = field(default_factory=list)
    ends: list = field(default_factory=list)
    energies: list = field(default_factory=list)
    events: list = field(default_factory=list)
    photon_ids: list = field(default_factory=list)

def transport_photons(pos, dirv, energy, geometry, rng,
                      grid=None, recorder=None) -> BatchResult:
```

記録ポイントはメインループ内に既にすべて揃っている:

- 飛行区間の始点 `o`・終点 `o + d*ds[:,None]`・エネルギー `e` は
  [transport.py:256-271](../vivemonte/transport.py) で計算済み。
  `recorder is not None` のときだけ、この`idx`（グローバル光子番号）と
  区間を追記する。
- 区間終端のイベント種別は同ループ内の分岐結果から決まる:
  `will_interact` かつ `is_photo/is_compt/is_rayl` → 各相互作用名、
  `~will_interact & escape` → "escape"、`~will_interact & ~escape` → "boundary"。
- **photon_id（元の光子番号 = `idx` の値）を必ず記録する**。
  可視化側で光子ごとの折れ線（ポリライン）に組み立てるのに必要。
- 実装はループ先頭で `if recorder is not None:` ブロックにまとめ、
  Noneのときのオーバーヘッドが分岐1回だけになるようにする。

`run_transport` には手を入れない（traceは独自にsample_source_photons →
transport_photonsを直接呼ぶ。__main__.pyのcmd_trace内で完結させる）。

### 実装 1-b: HTMLへの重ね描き（vivemonte/preview.py）

`scene_to_json` の返り値に `trajectories` キーを追加できるようにする:

```python
def scene_to_json(scene, trajectories=None) -> dict:
    ...
    return {"objects": ..., "beam": ..., "center": ..., "radius": ...,
            "warnings": ..., "trajectories": trajectories or []}
```

`trajectories` のJSON構造（光子ごとにまとめる。区間の羅列ではなく
ポリライン＋イベント列にするとJS側が単純になる）:

```json
[
  {"points": [[x,y,z], [x,y,z], ...],      // 始点＋各区間終点
   "energies": [120.0, 80.3, ...],          // 各区間のエネルギー[keV]
   "events": ["boundary", "compton", "photoelectric"]},  // 各区間終端
  ...
]
```

`_TEMPLATE` の拡張（既存の描画関数 `seg`/`project` をそのまま使う）:

- `draw()` 内、ビーム描画の後に軌跡描画ブロックを追加。
  `DATA.trajectories` が空なら何もしない（previewとの互換）。
- **色: エネルギーで連続変化**。`hue = 240 * (E / E_max)`（E_maxは
  全軌跡の最大エネルギー）とし、高エネルギー=青系→低エネルギー=赤系。
  `ctx.strokeStyle = `hsl(${hue},85%,60%)``。区間ごとに色を変える。
- **マーカー: イベント種別**。区間終端に:
  - photoelectric: 塗りつぶし円（半径3px）— 消滅点
  - compton: 白抜き円
  - rayleigh: 白抜き菱形
  - escape: 小さな「×」
  - boundary: マーカーなし（線が続くだけ）
- ヘッダーに `<label><input type="checkbox" id="ckTraj" checked> 軌跡</label>`
  を追加（trajectoriesが空のときはlabel自体を `display:none`）。
- 凡例にイベントマーカーの説明とエネルギー色スケールの説明を追加
  （trajectoriesが空のときは出さない）。

**注意**: `_TEMPLATE` は `.replace("__DATA__", data)` で埋め込むため、
JSON内に `</script>` 相当の文字列が入らないことは既存と同条件。
軌跡データは数値のみなので問題ない。

### 実装 1-c: CLI（vivemonte/__main__.py）

`cmd_trace(args)`:

1. `load_scene` → 検証エラーなら既存コマンドと同じ形式で終了
2. `sample_source_photons(src, n, rng)` → `transport_photons(..., recorder=rec)`
3. recorderの生データをphoton_idごとにグループ化してポリラインJSONに変換
   （このグループ化関数 `trajectories_to_json(recorder) -> list` は
   preview.pyかtransport.pyのどちらかに置き、**単体テスト対象にする**）
4. `write_html(scene, out, trajectories=...)`

### テスト（tests/test_trace.py 新規）

1. **区間の連続性**: 単一光子（n=1、seed固定）で、各ポリラインの
   `points[i+1]` が次区間の始点と一致する（軌跡が途切れていない）こと。
   注: 境界通過時の `pos += dirv * 1e-6` ナッジがあるため、
   許容誤差は `atol=1e-5` 程度にする。
2. **エネルギー単調性**: 任意の光子で `energies` は非増加
   （コンプトンでのみ減り、レイリー・境界では不変）。
3. **イベント整合性**: 最後のeventは必ず "photoelectric" か "escape"。
   途中に "photoelectric" が現れないこと。
4. **recorder=None の無影響**: 同一seedで recorder有無の
   `BatchResult` が完全一致（乱数消費が変わらないこと。
   recorderは乱数を一切引かない実装であることの担保）。
5. **HTMLスモーク**: cmd_trace相当を関数呼び出しで実行し、出力HTMLに
   `"trajectories"` と `ckTraj` が含まれ、previewで生成したHTMLには
   trajectoriesが空配列として入ることを確認。
6. **既存テスト無変更で全通過**（recorder追加が既存動作を変えない）。

---

## Deliverable 2: `vivemonte plot` — 線量マップの断面図

### CLI仕様

```bash
python -m vivemonte plot dose.npz [-o dose_maps.png] [--scene scene.yaml]
                                  [--quantity dose|h10] [--axis x|y|z] [--pos CM]
```

- `dose.npz` は既存 `run --dose-out` の出力（キー: `dose_per_history_Gy`,
  `h10_per_history_pSv`, 校正済みなら `dose_Gy`, `h10_pSv`,
  `origin_cm`, `voxel_size_cm`, `shape`）。
- `--quantity` 既定 `dose`。校正済みキー（`dose_Gy`）があればそれを優先し
  タイトルに「mAs校正済み [Gy]」、なければ per_history 値で「[Gy/history]」
  と単位を明記する。
- `--axis/--pos` 未指定時: **最大線量ボクセルを通る直交3断面**
  （axial/coronal/sagittal相当）を1枚のfigureに横並びで描く。
  指定時はその1断面のみ。
- `--scene` を渡すとジオメトリー輪郭をオーバーレイする（後述）。

### 実装（vivemonte/plotting.py 新規 + __main__.pyにcmd_plot）

- matplotlibは `matplotlib.use("Agg")`（既存cmd_xsと同じ流儀）。
- カラースケールは `LogNorm`。**vmin/vmaxの扱いが唯一の設計ポイント**:
  `vmax = data.max()`, `vmin = max(data[data>0].min(), vmax*1e-6)`。
  ゼロボクセル（線量なし）は `np.ma.masked_where(data<=0, ...)` で
  マスクし背景色（`cmap.set_bad`）にする。データが全ゼロなら
  エラーメッセージを出して終了コード1（クラッシュさせない）。
- 各断面のaxes: 物理座標軸[cm]（originとvoxel_sizeから`extent`を計算）、
  カラーバー付き、どの位置の断面か（例: `z = 140.0 cm`）をタイトルに。
- **ジオメトリー輪郭オーバーレイ**（`--scene`指定時のみ）:
  断面グリッドの各ピクセル中心で `Geometry.material_at` を評価し、
  材料名を整数IDに変換した2D配列に対して `plt.contour(levels=材料境界)`
  で輪郭線を引く。材料→整数IDは `np.unique(mat_2d, return_inverse=True)`
  で機械的に作る。輪郭は白（黒背景側でも見えるようlinewidth=0.7,
  alpha=0.8）。ラベルは不要（previewで確認済みの前提）。
- H*(10)は空気中に分布する量なのでマスクされにくく絵になる。
  dose(カーマ)は物体内に集中する。両方確認するのが本ワークフローの狙い。

### テスト（tests/test_plotting.py 新規）

1. 小さなnpz（例: 8×8×8、1ボクセルだけ非ゼロ）を組み立てて
   plot関数を呼び、pngファイルが生成され非ゼロサイズであること。
2. 全ゼロデータで例外にならず戻り値/終了コードでエラー通知すること。
3. 断面選択ロジック（最大ボクセルを通る3断面のインデックス計算）を
   純関数に切り出して直接テスト（`argmax` → unravel_indexの結果と一致）。
4. extent計算（origin/voxel_sizeから物理座標範囲）の数値照合。

---

## Deliverable 3: `/vive-check` スキル — 関門つき進行役

### 置き場所

**リポジトリ内** `.claude/skills/vive-check/SKILL.md`（このリポジトリは
公開GitHub管理なのでプロジェクト固有ワークフローとして版管理する。
KBの「ユーザーレベルに置く」ルールはKBの.claudeがgitignoreだから
であり、ここには当てはまらない）。`.gitignore` が `.claude/` を
除外していないことを実装時に確認すること。

### スキル本文の構成（実装を書かない。CLIを呼ぶ手順と判断だけ）

```markdown
---
name: vive-check
description: viveMonteのシーンを人間と確認しながら段階実行する
  （ジオメトリー確認→軌跡確認→本計算→結果確認の関門つきワークフロー）。
  「シーンを確認しながら実行して」「vive-checkして」などで起動。
---

# vive-check — 関門つき実行ワークフロー

引数: scene.yamlのパス（未指定なら examples/ 配下から確認する）

各関門で必ずユーザーの承認を待つ。承認前に次の段階を先行実行しない。

## 関門1: ジオメトリー確認
1. `.venv/bin/python -m vivemonte validate <scene>` — エラーなら修正提案して停止
2. `.venv/bin/python -m vivemonte preview <scene> -o <scene名>_preview.html`
3. `open <html>` でブラウザ表示し、AskUserQuestionで確認:
   「ジオメトリーは意図通りですか？」（はい / 修正が必要）
   修正指示があればscene.yamlを編集して関門1を再実行

## 関門2: ビーム・軌跡確認
1. `.venv/bin/python -m vivemonte trace <scene> -n 200 --seed 42 -o ..._trace.html`
2. `open <html>` で表示し、確認ポイントを言葉で添える:
   ビームの向き・照射野が意図通りか / 散乱が起きている場所は妥当か /
   遮蔽を光子が素通りしていないか
3. AskUserQuestionで承認を取る

## 関門3: 本計算
1. history数と解像度をユーザーに確認（既定: -n 1e6, --resolution 2）
   ※ 3e5 historiesで約15秒（線量グリッド併用時）。1e7超は時間を予告する
2. `.venv/bin/python -m vivemonte run <scene> -n <N> --seed <s> \
       --dose-grid --resolution <R> --dose-out <scene名>_dose.npz`
3. 実行結果のサマリ（吸収/脱出割合、mAs校正値など）を提示

## 関門4: 結果確認
1. `.venv/bin/python -m vivemonte plot <npz> --scene <scene> -o ..._maps.png`
2. 画像をReadで読んで内容を要約し、ユーザーに提示
3. 追加の断面や別の量（dose/h10）を求められたら--axis/--pos/--quantityで再描画
```

（上記はスキル本文の骨子。実装時はこの構成のままSKILL.mdに書く）

---

## 実装順序とコミット粒度

1. **コミット1**: TrajectoryRecorder + transport_photonsのrecorder引数
   ＋ trajectories_to_json ＋ tests/test_trace.py の1〜4
2. **コミット2**: preview.pyテンプレート拡張 + cmd_trace + テスト5
3. **コミット3**: plotting.py + cmd_plot + tests/test_plotting.py
4. **コミット4**: .claude/skills/vive-check/SKILL.md + README更新
   （使い方セクションにtrace/plotの例、実装済みリストに2項目追加）

各コミット前に `pytest tests/ -q` 全通過を確認。コミットメッセージは
日本語・本文で「何を・なぜ」（このリポジトリの慣習に従う）。

## 受け入れ基準

- [ ] 既存の全テストが**無変更で**通る（recorder追加の非破壊性）
- [ ] `vivemonte trace examples/chest_room.yaml -n 200` が10秒以内に
      完了し、ブラウザで軌跡・マーカー・凡例が表示される
- [ ] 同HTMLで軌跡チェックボックスをオフにすると既存previewと同等の表示
- [ ] `vivemonte run ... --dose-out d.npz` → `vivemonte plot d.npz --scene ...`
      で3断面＋輪郭＋カラーバーのpngが得られる
- [ ] recorder=None時のtransport_photonsが同一seedで従来と bit-identical
      な結果を返す（テスト4）
- [ ] スキル起動で関門1→4が順に進み、各関門でAskUserQuestionが出る

## 実装時の既知の罠（引き継ぎ）

- 境界通過時に `pos += dirv * 1e-6` のナッジがある（transport.py）。
  軌跡の連続性テストはこのぶんの許容誤差が要る。
- previewの `_TEMPLATE` はf-stringではなくraw文字列＋replace方式。
  `{}` を含むJSを追加しても壊れないが、`__DATA__` という文字列を
  JSコメント等に書かないこと。
- matplotlibのimportはコマンド関数内で行う（既存cmd_xsと同じ。
  起動時間とheadless環境への配慮）。
- `.npz` のキー有無で校正済みかを判定する（`"dose_Gy" in npz.files`）。
- 教訓ファイル（docs/lessons_learned.md）にある通り、物理量の絶対値を
  表示する箇所では単位と「per-historyか校正済みか」を必ず明記する。
