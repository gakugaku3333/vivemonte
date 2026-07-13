# viveMonte

診断X線領域（10〜150 keV）のモンテカルロ光子輸送計算システム。
AIがジオメトリー等を宣言的な `scene.yaml` で設定する前提で設計。

> ⚠️ 現状は教育・研究用。患者線量評価や遮蔽設計の実務判断に使う前に、
> PHITS等の確立コードとの相互検証を通すこと。

## セットアップ

```bash
python3 -m venv .venv
.venv/bin/pip install numpy pyyaml matplotlib xraylib pytest
```

## 使い方（すべて非対話CLI＝AIが叩ける）

```bash
# シーンを検証（物理サニティチェック付き）
.venv/bin/python -m vivemonte validate examples/chest_room.yaml

# ジオメトリーを3DプレビューHTMLに書き出し（外部依存なしの自己完結HTML）
.venv/bin/python -m vivemonte preview examples/chest_room.yaml -o preview.html

# 断面積カーブを描画
.venv/bin/python -m vivemonte xs water bone lead -o xs.png
```

## 断面積データの出どころ

| 量 | ソース | 用途 |
|---|---|---|
| μ/ρ・光電/コンプトン/レイリー内訳 | **xraylib**（EPDLベース、NIST XCOMと一致確認済み） | 輸送カーネルの自由行程・相互作用抽選 |
| μen/ρ（質量エネルギー吸収係数） | **NIST XAAMDI**（同梱CSV, Hubbell & Seltzer） | カーマ・吸収線量タリー |

`vivemonte/data/nist_xaamdi/` に20材料を同梱。再取得は
`.venv/bin/python scripts/fetch_nist_xaamdi.py`。
xraylibの `CS_Energy` はNIST公表μen/ρと最大約17%乖離するため使わない（`tests/`で検証）。

## テスト

```bash
.venv/bin/python -m pytest tests/ -q   # NIST公表値とのスポット照合
```

## 実装済み / 未実装

- [x] scene.yaml スキーマ＋検証器（自己修正ループ用の明確なエラー）
- [x] 3Dジオメトリープレビュー（自己完結HTML、回転/ズーム/断面）
- [x] 材料・断面積データ層（xraylib + NIST XAAMDI）
- [ ] 光子輸送カーネル（Woodcock tracking、バッチベクトル化）※次段
- [ ] ボクセル吸収線量・空間線量H*(10)タリー
- [ ] スペクトル生成（SpekPy統合）
- [ ] 解析解ベンチマーク（減弱・HVL・深部線量）
