---
name: vive-check
description: ChatCarloのシーンを人間と確認しながら段階実行する（ジオメトリー確認→軌跡確認→本計算→結果確認の関門つきワークフロー）。「シーンを確認しながら実行して」「vive-checkして」などで起動。
---

# vive-check — 関門つき実行ワークフロー

引数: scene.yamlのパス（未指定なら examples/ 配下から確認する）。

シーンがまだ存在しない・要件が曖昧な場合は、先に vive-interview スキルで
ヒアリングして scene.yaml を確定させてから、このスキルの関門1に入る。
vive-interviewから引き継いだ場合は `runs/<シーン名>/requirements.md`
（ヒアリングで確定した要件記録）のパスも受け取り、関門1のステージA監査に
渡す。

このスキルは**進行と判断だけ**を担う。実装（軌跡記録・断面図描画など）は
すべて `chatcarlo` CLI（`chatcarlo/transport.py` / `preview.py` / `plotting.py`、
テスト付き）にある。各関門で必ずユーザーの承認を待ち、承認前に次の段階を
先行実行しない。

**成果物の置き場所**: preview/trace/dose/mapsの出力ファイルはすべて
`runs/<シーン名>/`（`.gitignore`済みの作業用ディレクトリ）にまとめる。
リポジトリ直下に散らかさない。関門1の最初に `mkdir -p runs/<シーン名>` する
（vive-interviewから引き継いだ場合はすでに存在するはず）。最終報告に残す
価値のある成果物（学会抄録用の実演結果など）は、`docs/egs5_crosscheck/<name>/`
のように意図的にコピーして追跡する（`runs/`はgitignore対象で自動では残らない）。

## 関門1: ジオメトリー確認

1. `mkdir -p runs/<シーン名>`（未作成なら）。
   `.venv/bin/python -m chatcarlo validate <scene>` — エラーが出たら
   scene.yamlの修正案を提示して停止する（先へ進まない）。
2. `.venv/bin/python -m chatcarlo preview <scene> -o runs/<シーン名>/preview.html`
3. vive-auditスキルでステージA（シーン監査）を実施し、監査報告を提示する
   （ユーザーが監査を不要と言った場合は省略可）。`runs/<シーン名>/requirements.md`
   があれば必ずそのパスも監査依頼に含める（ユーザー意図との一致検証に使う）。
   差し戻し所見があればscene.yamlを修正して再監査してから次へ。
4. `open <html>` でブラウザ表示し、AskUserQuestionで確認する:
   「ジオメトリーは意図通りですか？」（はい / 修正が必要）
   修正指示があればscene.yamlを編集し、関門1を最初からやり直す。

## 関門2: ビーム・軌跡確認

1. `.venv/bin/python -m chatcarlo trace <scene> -n 200 --seed 42 -o runs/<シーン名>/trace.html`
   （n=200は既定。2000を超えると自動的にクランプされ警告が出る。10秒以内に
   終わるはず）
2. `open <html>` で表示する。確認ポイントを言葉で添える:
   - ビームの向き・照射野が意図通りか
   - 散乱（○コンプトン・◇レイリー）が起きている場所は妥当か
   - 遮蔽（lead等）を光子が素通りしていないか（●光電吸収が遮蔽内で
     十分起きているか）
3. AskUserQuestionで承認を取る（はい / 修正が必要）。
   修正指示があればscene.yamlを編集し、関門1からやり直す。

## 関門3: 本計算

1. history数と解像度をユーザーに確認する（既定案: `-n 1e6 --resolution 2`）。
   目安: 3e5 historiesで線量グリッド併用時に約15秒。1e7を超える場合は
   実行前に「時間がかかります」と予告する。
2. `.venv/bin/python -m chatcarlo run <scene> -n <N> --seed <s> \
       --dose-grid --resolution <R> --dose-out runs/<シーン名>/dose.npz`
3. vive-auditスキルでステージB（実行結果監査）を実施する。run の標準出力
   全文を監査官に渡し、独立オーダー照合（Beer-Lambert / SpekPyカーマ）まで
   済ませてもらう。
4. 実行結果のサマリ（吸収/脱出割合、mAs校正値など、コマンドの標準出力）を
   監査報告書とともにユーザーに提示する。

## 関門4: 結果確認

1. `.venv/bin/python -m chatcarlo plot runs/<シーン名>/dose.npz --scene <scene> \
       -o runs/<シーン名>/maps.png`
   （既定は最大値ボクセルを通る3断面。`--axis x|y|z --pos <cm>` で特定断面、
   `--quantity h10` でH*(10)マップに切り替えられる）
2. 画像をReadツールで読んで内容を要約し、ユーザーに提示する。
3. 追加の断面や別の量を求められたら `--axis`/`--pos`/`--quantity` を変えて
   再描画する。
4. 数値を含む最終報告文を書いたら、送信前にvive-auditスキルで
   ステージC（報告文監査）にかける。

## 単位・校正済みかどうかの明記（重要）

線量・H*(10)の絶対値を口頭で報告するときは、必ず「1 historyあたり」か
「校正済みの実測相当値」（mAs基準またはCTDIvol基準）かを明示すること。
`chatcarlo plot` のカラーバーラベル（`Gy/history` か `Gy (calibrated)` 等）にも
既にこの区別が入っている。過去に絶対値校正の桁を間違えたことがある
（`docs/lessons_learned.md` 参照）ため、数値だけを鵜呑みにせず単位表記を
確認してから報告する。
