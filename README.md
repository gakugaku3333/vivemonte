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

# 光子輸送を実行（材料別吸収エネルギー・吸収/脱出割合を表示）
.venv/bin/python -m vivemonte run examples/chest_room.yaml -n 1e6 --seed 42
```

## 輸送カーネルの設計

ボクセル化はしていない。box/cylinder/sphere の解析プリミティブのまま、
光子ごとに「次の材料境界までの距離」を都度計算して進む**解析面トラッキング**
方式（`vivemonte/geometry.py` + `vivemonte/transport.py`）。

- 部屋スケール（数百cm）× 薄い遮蔽（数cm）を素直に表現でき、ボクセル解像度と
  メモリのトレードオフが発生しない
- 各区間内は材料が均質なので μ が一定 → Woodcock delta-trackingの仮想衝突が
  不要（空気の広い空間で無駄な計算をしない）
- 重なりはリスト順で後勝ち。物体を跨がない開放空間は background（既定 air）
- 相互作用は光電/コンプトン（Klein-Nishina, Kahn型棄却法）/レイリー
  （簡易Thomson近似）の3種。電子飛程を無視するカーマ近似で局所吸収を計算
- `tests/test_transport.py` で単色鉛筆ビームの一次透過率をBeer-Lambert則
  （exp(-μt)）と直接照合

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
- [x] 光子輸送カーネル（解析面トラッキング、numpyバッチベクトル化）
- [x] 解析解ベンチマーク（一次透過率とBeer-Lambert則の照合）
- [ ] ボクセル吸収線量・空間線量H*(10)タリー（グリッドは輸送とは独立に敷く）※次段
- [ ] スペクトル生成（SpekPy統合。現状はKramers則＋Al濾過減弱の簡易近似）
- [ ] レイリー散乱の原子形状因子（現状はThomson型 (1+cos²θ)/2 の粗い近似）
- [ ] 光電吸収時の蛍光X線を無視（局所吸収として扱っている）
