# viveMonte

診断X線領域（10〜150 keV）のモンテカルロ光子輸送計算システム。
AIがジオメトリー等を宣言的な `scene.yaml` で設定する前提で設計。

> ⚠️ 現状は教育・研究用。患者線量評価や遮蔽設計の実務判断に使う前に、
> PHITS等の確立コードとの相互検証を通すこと。

## セットアップ

```bash
python3 -m venv .venv
.venv/bin/pip install numpy pyyaml matplotlib xraylib pytest spekpy
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

# ボクセル吸収線量タリーも計算し.npzに書き出す
.venv/bin/python -m vivemonte run examples/chest_room.yaml -n 1e6 --seed 42 \
    --dose-grid --resolution 5 --dose-out dose.npz
```

scene.yamlの`source.mas`（管電流時間積 [mAs]）を指定すると、SpekPyの絶対
フルエンス計算で実際に照射野を通過する光子数を求め、`Gy/history`・
`pSv/history`の相対値を実際の`Gy`・`pSv`（1回の撮影あたり）に校正する。
未指定時は相対値（history単位）のみを出力する。

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
  （原子形状因子F(Z,q)込み、xraylib.FF_Rayl）の3種。電子飛程を無視する
  カーマ近似で局所吸収を計算
- `tests/test_transport.py` で単色鉛筆ビームの一次透過率をBeer-Lambert則
  （exp(-μt)）と直接照合

ボクセル吸収線量・H*(10)（`vivemonte/tally.py`）はタリー専用の一様グリッドとして
輸送ジオメトリーとは独立に敷く。相互作用点ごとに集計する
collision estimator（`energy_deposited`）と、飛程積分によるカーマ
track-length estimator（グリッド）の**2つの独立推定量が統計誤差内で一致する**
ことをテストで相互検証している（`tests/test_tally.py`）。H*(10)はフルエンス
ベースの防護量でカーマとは異なる量のため、飛程積分をボクセル体積で正規化して
求める（`VoxelGrid.h10_map_pSv`）。線量は既定で`Gy/history`・`pSv/history`
単位で出力し、`source.mas`指定時はSpekPyの絶対フルエンスで実際の
`Gy`・`pSv`に校正する（`vivemonte/transport.py: photon_count_through_field`）。

## 断面積データの出どころ

| 量 | ソース | 用途 |
|---|---|---|
| μ/ρ・光電/コンプトン/レイリー内訳 | **xraylib**（EPDLベース、NIST XCOMと一致確認済み） | 輸送カーネルの自由行程・相互作用抽選 |
| μen/ρ（質量エネルギー吸収係数） | **NIST XAAMDI**（同梱CSV, Hubbell & Seltzer） | カーマ・吸収線量タリー |
| h*(10)/Φ（周辺線量当量換算係数） | **ICRP Publication 74 / ICRU Report 57**（同梱CSV） | 空間線量H*(10)タリー |

`vivemonte/data/nist_xaamdi/` に20材料を同梱。再取得は
`.venv/bin/python scripts/fetch_nist_xaamdi.py`。
xraylibの `CS_Energy` はNIST公表μen/ρと最大約17%乖離するため使わない（`tests/`で検証）。

`vivemonte/data/h_star_10/photons_icrp74.csv` にICRP74公表のh*(10)/Φ表を同梱。
再取得は `.venv/bin/python scripts/fetch_h_star_10.py`（OpenMCプロジェクト同梱の
転記データを取得、値自体はICRP74）。診断X線領域ではH*(10)は個人線量当量
Hp(10,0°)・実効線量E(AP)と数値的にほぼ一致することが知られている
（Otto, JINST 2019, arXiv:1906.05411 — Eph ≤ 6 MeVで成立）。

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
- [x] ボクセル吸収線量タリー（グリッドは輸送とは独立に敷く。track-length estimator）
- [x] 空間線量H*(10)タリー（ICRP74公表のh*(10)/Φ表、フルエンスのtrack-length estimator）
- [x] 線源のフォトン数校正（`source.mas`指定時、SpekPyの絶対フルエンスで
      Gy/history・pSv/historyをGy・pSvに換算。フィールド内フルエンスは
      中心軸上の値で一様と近似する一次近似）
- [x] スペクトル生成（SpekPy統合。タングステン陽極、既定陽極角12度、
      `source.anode_angle_deg`で変更可。SpekPy未インストール時はKramers則＋
      Al濾過減弱の簡易近似にフォールバック）
- [x] レイリー散乱の原子形状因子（xraylib.FF_Rayl、EPDLベース。化合物は
      質量分率×元素別断面積で構成元素を抽選してから角度分布を決める）
- [ ] 光電吸収時の蛍光X線を無視（局所吸収として扱っている）
