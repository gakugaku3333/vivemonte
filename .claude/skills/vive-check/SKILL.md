---
name: vive-check
description: viveMonteのシーンを人間と確認しながら段階実行する（ジオメトリー確認→軌跡確認→本計算→結果確認の関門つきワークフロー）。「シーンを確認しながら実行して」「vive-checkして」などで起動。
---

# vive-check — 関門つき実行ワークフロー

引数: scene.yamlのパス（未指定なら examples/ 配下から確認する）。

シーンがまだ存在しない・要件が曖昧な場合は、先に vive-interview スキルで
ヒアリングして scene.yaml を確定させてから、このスキルの関門1に入る。

このスキルは**進行と判断だけ**を担う。実装（軌跡記録・断面図描画など）は
すべて `vivemonte` CLI（`vivemonte/transport.py` / `preview.py` / `plotting.py`、
テスト付き）にある。各関門で必ずユーザーの承認を待ち、承認前に次の段階を
先行実行しない。

## 関門1: ジオメトリー確認

1. `.venv/bin/python -m vivemonte validate <scene>` — エラーが出たら
   scene.yamlの修正案を提示して停止する（先へ進まない）。
2. `.venv/bin/python -m vivemonte preview <scene> -o <scene名>_preview.html`
3. `open <html>` でブラウザ表示し、AskUserQuestionで確認する:
   「ジオメトリーは意図通りですか？」（はい / 修正が必要）
   修正指示があればscene.yamlを編集し、関門1を最初からやり直す。

## 関門2: ビーム・軌跡確認

1. `.venv/bin/python -m vivemonte trace <scene> -n 200 --seed 42 -o <scene名>_trace.html`
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
2. `.venv/bin/python -m vivemonte run <scene> -n <N> --seed <s> \
       --dose-grid --resolution <R> --dose-out <scene名>_dose.npz`
3. 実行結果のサマリ（吸収/脱出割合、mAs校正値など、コマンドの標準出力）を
   ユーザーに提示する。

## 関門4: 結果確認

1. `.venv/bin/python -m vivemonte plot <scene名>_dose.npz --scene <scene> \
       -o <scene名>_maps.png`
   （既定は最大値ボクセルを通る3断面。`--axis x|y|z --pos <cm>` で特定断面、
   `--quantity h10` でH*(10)マップに切り替えられる）
2. 画像をReadツールで読んで内容を要約し、ユーザーに提示する。
3. 追加の断面や別の量を求められたら `--axis`/`--pos`/`--quantity` を変えて
   再描画する。

## 単位・校正済みかどうかの明記（重要）

線量・H*(10)の絶対値を口頭で報告するときは、必ず「1 historyあたり」か
「mAs校正済みの実測相当値」かを明示すること。`vivemonte plot` のカラーバー
ラベル（`Gy/history` か `Gy (mAs-calibrated)` 等）にも既にこの区別が
入っている。過去に絶対値校正の桁を間違えたことがある
（`docs/lessons_learned.md` 参照）ため、数値だけを鵜呑みにせず単位表記を
確認してから報告する。
