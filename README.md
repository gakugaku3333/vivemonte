# ChatCarlo

診断X線領域（10〜150 keV）のモンテカルロ光子輸送計算システム。
AIがジオメトリー等を宣言的な `scene.yaml` で設定する前提で設計。

> ⚠️ 現状は教育・研究用。患者線量評価や遮蔽設計の実務判断に使う前に、
> EGS5等の確立コードとの相互検証を通すこと。

## セットアップ

```bash
python3 -m venv .venv
.venv/bin/pip install numpy pyyaml matplotlib xraylib pytest spekpy scipy
```

## 使い方（すべて非対話CLI＝AIが叩ける）

```bash
# シーンを検証（物理サニティチェック付き）
.venv/bin/python -m chatcarlo validate examples/chest_room.yaml

# ジオメトリーを3DプレビューHTMLに書き出し（外部依存なしの自己完結HTML）
.venv/bin/python -m chatcarlo preview examples/chest_room.yaml -o preview.html

# 断面積カーブを描画
.venv/bin/python -m chatcarlo xs water bone lead -o xs.png

# 光子輸送を実行（材料別吸収エネルギー・吸収/脱出割合を表示）
.venv/bin/python -m chatcarlo run examples/chest_room.yaml -n 1e6 --seed 42

# ボクセル吸収線量タリーも計算し.npzに書き出す
.venv/bin/python -m chatcarlo run examples/chest_room.yaml -n 1e6 --seed 42 \
    --dose-grid --resolution 5 --dose-out dose.npz

# 光子軌跡を3D可視化（小history。既存previewのHTMLに軌跡を重ね描き）
.venv/bin/python -m chatcarlo trace examples/chest_room.yaml -n 200 --seed 42 -o trace.html

# 線量/H*(10)マップの断面図（既定: 最大値ボクセルを通る3断面＋ジオメトリー輪郭）
.venv/bin/python -m chatcarlo plot dose.npz --scene examples/chest_room.yaml -o maps.png
```

`.claude/skills/vive-check/` に、上記4コマンドを「ジオメトリー確認→軌跡確認→
本計算→結果確認」の順で人間の承認を挟みながら進める進行役スキルを同梱している
（「シーンを確認しながら実行して」等で起動）。

`.claude/skills/vive-interview/` に、シーン未作成の段階で「どんなシミュレーションを
したいのか」を目的→撮影条件→ジオメトリー→計算条件の順にヒアリングし、不明瞭な
部分を質問で詰めてからscene.yaml草案を作る聞き取りスキルを同梱している
（「〜の被ばくを計算したい」等の相談で起動、完成後はvive-checkに引き継ぐ）。

`examples/` に検査種別ごとのベースシーンを同梱している（ヒアリング時の雛形）:

| ファイル | 検査種別 |
|---|---|
| `chest_room.yaml` | 一般撮影・立位胸部（撮影室の空間線量・遮蔽） |
| `abdomen_table.yaml` | 一般撮影・臥位腹部（テーブルBucky） |
| `ct_room.yaml` | CT（連続ガントリー回転＋ヘリカルの位相平均近似。⚠️ボウタイ・mA変調未対応） |
| `angio_room.yaml` | アンギオ/IVR・透視（術者被ばく: 防護衣＋天吊り遮蔽つき） |
| `portable_ward.yaml` | ポータブル撮影（隣ベッド患者への散乱線） |

scene.yamlの`source.mas`（管電流時間積 [mAs]）を指定すると、SpekPyの絶対
フルエンス計算で実際に照射野を通過する光子数を求め、`Gy/history`・
`pSv/history`の相対値を実際の`Gy`・`pSv`（1回の撮影あたり）に校正する。
未指定時は相対値（history単位）のみを出力する。

> ⚠️ `--dose-grid`が報告する「最大吸収線量」「最大H*(10)」は、多くの場合
> 背景（空気）ボクセルに落ちる（線源近傍の1/r²発散、材料境界での後方散乱等）。
> その場合は`[警告]`が表示される。実在する位置（患者表面・操作者位置等）の
> 被ばく評価には、その位置に直接細かいグリッドを敷いて計算すること
> （詳細は[docs/lessons_learned.md](docs/lessons_learned.md)）。

## 輸送カーネルの設計

ボクセル化はしていない。box/cylinder/sphere の解析プリミティブのまま、
光子ごとに「次の材料境界までの距離」を都度計算して進む**解析面トラッキング**
方式（`chatcarlo/geometry.py` + `chatcarlo/transport.py`）。

- 部屋スケール（数百cm）× 薄い遮蔽（数cm）を素直に表現でき、ボクセル解像度と
  メモリのトレードオフが発生しない
- 各区間内は材料が均質なので μ が一定 → Woodcock delta-trackingの仮想衝突が
  不要（空気の広い空間で無駄な計算をしない）
- 重なりはリスト順で後勝ち。物体を跨がない開放空間は background（既定 air）
- 相互作用は光電/コンプトン（Klein-Nishina, Kahn型棄却法）/レイリー
  （原子形状因子F(Z,q)込み、xraylib.FF_Rayl）の3種。電子飛程を無視する
  カーマ近似で局所吸収を計算。光電吸収はK殻蛍光X線の放出を抽選し
  （`physics.fluorescence`、既定true）、放出時はエネルギーの一部を
  持ち出す光子として輸送を続ける（詳細は下記「実装済み」参照）
- `tests/test_transport.py` で単色鉛筆ビームの一次透過率をBeer-Lambert則
  （exp(-μt)）と直接照合

ボクセル吸収線量・H*(10)（`chatcarlo/tally.py`）はタリー専用の一様グリッドとして
輸送ジオメトリーとは独立に敷く。相互作用点ごとに集計する
collision estimator（`energy_deposited`）と、飛程積分によるカーマ
track-length estimator（グリッド）の**2つの独立推定量が統計誤差内で一致する**
ことをテストで相互検証している（`tests/test_tally.py`）。H*(10)はフルエンス
ベースの防護量でカーマとは異なる量のため、飛程積分をボクセル体積で正規化して
求める（`VoxelGrid.h10_map_pSv`）。線量は既定で`Gy/history`・`pSv/history`
単位で出力し、`source.mas`指定時はSpekPyの絶対フルエンスで実際の
`Gy`・`pSv`に校正する（`chatcarlo/source.py: photon_count_through_field`）。

## 断面積データの出どころ

| 量 | ソース | 用途 |
|---|---|---|
| μ/ρ・光電/コンプトン/レイリー内訳 | **xraylib**（EPDLベース、NIST XCOMと一致確認済み） | 輸送カーネルの自由行程・相互作用抽選 |
| μen/ρ（質量エネルギー吸収係数） | **NIST XAAMDI**（同梱CSV, Hubbell & Seltzer） | カーマ・吸収線量タリー |
| h*(10)/Φ（周辺線量当量換算係数） | **ICRP Publication 74 / ICRU Report 57**（同梱CSV） | 空間線量H*(10)タリー |

`chatcarlo/data/nist_xaamdi/` に20材料を同梱。再取得は
`.venv/bin/python scripts/fetch_nist_xaamdi.py`。
xraylibの `CS_Energy` はNIST公表μen/ρと最大約17%乖離するため使わない（`tests/`で検証）。

`chatcarlo/data/h_star_10/photons_icrp74.csv` にICRP74公表のh*(10)/Φ表を同梱。
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
- [x] CTガントリー回転＋ヘリカルの位相平均近似（`source.rotation`。角度は
      連続一様抽選、`scan_length_cm`で体軸方向にも一様分布。`n_angles`指定で
      離散投影モダリティも模擬可）
- [x] cone照射野（`field.shape: cone` + `diameter_cm`。円錐内で方向を立体角
      一様に抽選、CBCT・広角視野向け。rect照射野のSID面一様抽選は視野端で
      立体角あたりcos³θ分過剰になるが、coneは実管球の等方放出に対応し
      この偏りがない。回転・ヘリカルとも合成可）
- [x] ヒール効果（`source.heel_effect: true` + `anode_direction`。SpekPyの
      軸外スペクトルをヒール軸に沿ってビン化し、陽極側の強度低下を棄却
      サンプリング、線質硬化をビン別エネルギー抽選で再現。陽極カットオフを
      超える照射野は警告。rect/cone・回転・ヘリカルと合成可、mAs校正は
      照射野の面平均フルエンスに切替）
- [x] CTDIvol基準の絶対線量校正（`source.ctdi_vol_mGy`。コンソール表示値を
      アンカーに、標準CTDIファントム（PMMA Ø32/16cm）の内部シミュレーションで
      実効光子数を求める。ボウタイ・mA変調・ピッチが実測値に折り込まれるため
      CTではmAs校正より汎用的。`chatcarlo/ctdi.py`、mAs経路との桁の相互検証
      テスト付き）
- [x] スペクトル生成（SpekPy統合。タングステン陽極、既定陽極角12度、
      `source.anode_angle_deg`で変更可。SpekPy未インストール時はKramers則＋
      Al濾過減弱の簡易近似にフォールバック）
- [x] レイリー散乱の原子形状因子（xraylib.FF_Rayl、EPDLベース。化合物は
      質量分率×元素別断面積で構成元素を抽選してから角度分布を決める）
- [x] 光子軌跡の3D可視化（`chatcarlo trace`。既存previewテンプレートを拡張し
      軌跡なしでは従来のpreviewと同一表示になる設計。乱数を消費しない
      opt-inレコーダーで既存輸送結果に非破壊）
- [x] 線量/H*(10)マップの断面図（`chatcarlo plot`。最大値ボクセルを通る
      3断面が既定、`--scene`でジオメトリー輪郭を重ねられる）
- [x] 人間と確認しながら進める進行役スキル（`.claude/skills/vive-check/`。
      ジオメトリー確認→軌跡確認→本計算→結果確認の4関門、各関門でユーザー承認）
- [x] EGS5相互検証パイプライン（`.claude/skills/vive-crosscheck/`。計画・経緯:
      [docs/plan_egs5_crosscheck.md](docs/plan_egs5_crosscheck.md)。Phase 1
      〈一次透過率〉は合格、Phase 2a〈BSF〉は保留（自由空気カーマ測定の食い違いが
      未解消）、Phase 2b〈PDD/側方プロファイル〉は事前登録基準に対して不合格
      （EGS5−ChatCarlo系統差−1.71%、原因の大半を特定済み・調査継続中）。
      PHITS導入は取りやめ、EGS5を継続使用する）
- [x] 光電吸収時のK殻蛍光X線（`physics.fluorescence`、既定true。L殻・カスケードは
      非対応。線エネルギー5 keV未満はK端超でも局所吸収扱い。データは全てxraylib
      〈EdgeEnergy/FluorYield/RadRate/LineEnergy/CS_Photo_Partial、EPDLベース〉。
      計画: [docs/plan_fluorescence.md](docs/plan_fluorescence.md)）
