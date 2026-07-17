# 計画: 蛍光X線（特性X線）の実装

作成日: 2026-07-17 / ステータス: 未着手 / 実行担当: Claude (Sonnet可)

この計画書は、設計判断をすべて事前に固定してあり、実行者は各Phaseのタスクを
上から順に実施するだけでよい。判断に迷う箇所が出たら**勝手に設計変更せず**、
本文書の「設計判断（確定事項）」に立ち返るか、ユーザーに確認すること。

## 背景と目的

現状の光電吸収は「光子消滅・エネルギー全量をその場で局所吸収」
（[transport.py](../chatcarlo/transport.py) の `is_photo` 分岐）。実際には、
光電吸収で生じた内殻空孔が埋まる際に蛍光X線（特性X線）が放出され、
そのエネルギーは相互作用点から**持ち出される**。

診断X線領域（10–150 keV）でこれが効くのは高Z材料のK殻蛍光:

| 元素 | K吸収端 | Kα1 | K蛍光収率 ω_K | 影響 |
|---|---|---|---|---|
| Pb (Z=82) | 88.0 keV | 75.0 keV | ~0.96 | 遮蔽計算。88 keV超の光子が鉛に吸収されると約96%の確率で75–85 keVの透過力の高い蛍光X線が再放出される |
| W (Z=74) | 69.5 keV | 59.3 keV | ~0.95 | ターゲット・コリメータ |
| Cu/Fe/Ca | 7–9 keV | 6.4–8.0 keV | ~0.3–0.4 | 蛍光光子の自material内での平均自由行程がμmオーダーで実質局所吸収（後述のカットオフで処理） |

つまり**鉛遮蔽の透過線量・遮蔽背後のH*(10)を現状は過小評価しうる**
（K蛍光の持ち出し分を全部その場に落としているため、遮蔽体内で線量過大・
背後で過小）。EGS5相互検証（IEDGFL=1条件）でもこの差が系統差として出る。

副次効果: 現状は collision estimator（光電で全量計上）と track-length kerma
estimator（NIST μen/ρベース。μen/ρは平均蛍光放出分を**既に差し引いた**係数）
が高Z材料で原理的に食い違っている。蛍光実装で collision 側も蛍光分を差し引く
ようになり、両推定量の整合が改善するはず（Phase 3で確認する）。

## 設計判断（確定事項）

1. **K殻蛍光のみ実装する。L殻はやらない。**
   根拠: Pb のL蛍光（10.5–12.6 keV）は鉛中の平均自由行程が数μmで実質逃げない。
   診断領域で輸送に効くのはK蛍光だけ。カスケード（K空孔→L蛍光の連鎖）も
   同じ理由で無視し、残余エネルギーは局所吸収する。
2. **二次粒子スタックは作らない。その場で「光子を蛍光光子に置き換える」。**
   光電吸収1回が生む蛍光光子は高々1個なので、配列スロットを再利用できる:
   エネルギーを蛍光線エネルギーに書き換え、方向を等方に再抽選し、τを引き直して
   輸送続行。`E - E_line` を局所吸収に計上する。既存のベクトル化構造
   （aliveマスク方式）を一切変えずに済む — これが本設計の核。
3. **蛍光を放出しない場合**（非K殻イオン化 / K殻でも非放射遷移(1-ω_K) /
   線エネルギーがカットオフ未満）は現状どおり全量局所吸収・光子消滅。
4. **低エネルギーカットオフ**: 蛍光線エネルギーが **5 keV 未満**なら生成せず
   局所吸収（Ca/Fe/Cu等のK線が該当。自material内mfpがμmオーダーのため）。
   定数 `_FLUOR_CUTOFF_KEV = 5.0` を physics.py に置く。
   既存輸送の（コンプトンで下がった）低エネルギー光子の扱いは**変えない**。
5. **データはすべて xraylib から取る**（EPDLベース。手打ち転記禁止 —
   CLAUDE.md のデータ来歴ルール参照）:
   - K吸収端: `xraylib.EdgeEnergy(Z, xraylib.K_SHELL)` [keV]
   - K殻イオン化確率: `xraylib.CS_Photo_Partial(Z, xraylib.K_SHELL, E) / xraylib.CS_Photo(Z, E)`
   - K蛍光収率: `xraylib.FluorYield(Z, xraylib.K_SHELL)`
   - 線の分岐比と線エネルギー: `xraylib.RadRate(Z, line)` / `xraylib.LineEnergy(Z, line)`
     で line ∈ {`xraylib.KL2_LINE`(Kα2), `xraylib.KL3_LINE`(Kα1),
     `xraylib.KM2_LINE`, `xraylib.KM3_LINE`(Kβ1)} の4線。RadRateの和で規格化する
     （4線でK放射遷移の99%以上を占める。より弱い線は無視）。
   - 注意: `RadRate`/`LineEnergy` は該当線が存在しないZでエラーまたは0を返す
     ことがある。try/except で0扱いにし、有効線だけで規格化すること。
6. **蛍光光子の角度分布は等方**（内殻空孔の緩和は入射方向とほぼ無相関。
   EGS5等の標準的扱いと同じ）。
7. **scene.yaml にトグルを追加**: `physics.fluorescence: true|false`、
   **デフォルト true**（物理的により正しい側を既定にする）。過去のEGS5
   相互検証（IEDGFL未指定=蛍光なしで実施した分）を再現したいときに false にする。
   `physics` キー自体が新設なので scene.py のスキーマ検証に追加が必要。
8. **乱数の再現性**: 蛍光判定・線抽選は既存の輸送用 `rng` をそのまま使う
   （spawnしない）。同一seedの結果は蛍光実装前と**変わってよい**が、
   実装後の同一seed再現性（`tests/test_reproducibility.py`）は維持すること。

## 非目標（やらないこと）

- L殻蛍光・カスケード緩和・Auger電子（電子輸送自体がカーマ近似で無い）
- コンプトン後のドップラー広がり・電離緩和（束縛コンプトンのS(Z,q)までが現状）
- 制動放射
- xraylib以外のデータソース追加

## フェーズ計画

### Phase 1 — データ層（materials.py）

**タスク**:
1. `materials.py` に以下を追加:
   - `photo_element_weights(material, energies_keV)`:
     既存の `compton_element_weights` / `rayleigh_element_weights` と完全同型。
     断面積関数だけ `xraylib.CS_Photo` に変える。docstringに「光電相互作用が
     どの構成元素で起きたかの抽選重み（蛍光X線サンプリング用）」と書く。
   - `fluorescence_k_data(z: int)`（`@functools.lru_cache`）:
     戻り値 `(edge_keV, omega_k, line_energies: np.ndarray, line_probs: np.ndarray)`。
     line_probs は RadRate を有効線だけで規格化した累積でない生の確率
     （和=1）。有効線が0本（軽元素）なら `omega_k=0.0` 扱いで返す。
2. `tests/test_fluorescence.py` を新規作成し、データ層のスポットチェック:
   - Pb: K端 88.00±0.1 keV、Kα1 74.97±0.1 keV、ω_K 0.96±0.02
   - W: K端 69.53±0.1 keV、Kα1 59.32±0.1 keV
   - 各Zで line_probs の和が1（有効線がある場合）
   - `photo_element_weights("water", [50.0])` の重みの和が1、酸素が支配的
     （>0.9。水の光電は実質酸素）

**完了条件**: 上記テストが通る。既存テスト全体も無傷
（`.venv/bin/python -m pytest tests/ -q`）。

### Phase 2 — サンプリング層（physics.py）と輸送組み込み（transport.py）

**タスク**:
1. `physics.py` に追加:
   ```
   _FLUOR_CUTOFF_KEV = 5.0

   def sample_fluorescence(materials, e_keV, rng):
       """光電吸収イベント群に対し蛍光放出を抽選する。

       戻り値: (emit: bool配列, e_line: float配列)
       emit=True の光子は e_line [keV] のK蛍光光子を等方放出する。
       手順（各光子ごと、ベクトル化はmaterial_groups/Zグループ単位でよい）:
       1. photo_element_weights で吸収元素Zを抽選（sample_compton_elementと同型）
       2. fluorescence_k_data(Z) を引く。E <= K端 なら emit=False
       3. K殻イオン化確率 CS_Photo_Partial/CS_Photo（Eごとに評価）で棄却
       4. ω_K で棄却
       5. 線を line_probs で抽選 → e_line。e_line < _FLUOR_CUTOFF_KEV なら emit=False
       """
   ```
   実装メモ: `CS_Photo_Partial` はスカラー関数なので、Zグループ内でも
   エネルギーごとにループが要る。光電イベントは1バッチ数千件程度なので
   Pythonループで許容（既存のレイリーF(Z,q)ループと同じ割り切り）。
   K端直下のEで `CS_Photo_Partial` を呼ぶとエラーになるZがある —
   E <= K端 の判定を**先に**行い、その光子では呼ばないこと。
2. `transport.py` の `is_photo` 分岐を変更。現在の
   「全員 `_deposit` → alive=False」を以下に置換:
   - scene由来のフラグ `fluorescence_enabled`（`transport_photons` の引数に追加、
     デフォルト True。`run_transport` が scene の `physics.fluorescence` を渡す）
     が False なら現行動作のまま。
   - True なら `sample_fluorescence(mat_i[is_photo], e_i[is_photo], rng)` を呼び、
     - emit=False の光子: 現行どおり全量 `_deposit`、absorbed=True
     - emit=True の光子: `_deposit` は `e - e_line` のみ。
       `energy[photo_idx_emit] = e_line`、
       `dirv[photo_idx_emit] =` 等方抽選（`scatter_direction` は使わず、
       cosθ〜U(-1,1)・φ〜U(0,2π) から直接単位ベクトルを作る小関数
       `isotropic_direction(n, rng)` を physics.py に追加して使う）、
       `tau[photo_idx_emit] = -np.log(rng.random(...))`、
       alive のまま、`n_scatter += 1`。absorbed にはしない。
   - recorder には emit=True の光子のイベント名 `"fluorescence"` を記録
     （trace.html の色分けに新イベント種が増える。trajectory.py/テンプレートが
     未知イベント名で落ちないか確認し、必要なら凡例に1色追加）。
3. `scene.py`: トップレベル `physics:` キー（省略可、dict）を許可し、
   `fluorescence`（bool、デフォルトtrue）だけを検証。未知キーは既存の
   スタイルに合わせてエラーにする。`run_transport` 側で
   `scene.raw.get("physics", {}).get("fluorescence", True)` を参照。
4. テスト追加（`tests/test_fluorescence.py` に追記）:
   - **エネルギー保存**: 100 keV単色平行ビーム→鉛スラブ（十分厚: 2mm以上）の
     小シーンを `transport_photons` 直叩きで実行し、
     `sum(energy_deposited) + sum(escaped光子の final_energy)` が
     入射総エネルギーと一致（相対1e-9、蛍光ON/OFF両方で成立すること）
   - **蛍光の発生率**: 100 keV→薄い鉛箔で、K端超なので蛍光イベントが
     発生していること（recorder経由で `"fluorescence"` イベント数>0）。
     87 keV（K端未満）では蛍光イベント数=0
   - **脱出光子のスペクトル**: 100 keV→鉛薄箔からの脱出光子のうち、
     72–88 keV帯に蛍光ピーク由来の光子が存在する
     （final_energy のヒストグラムで該当ビンが非蛍光時より増える、程度の緩い判定）
   - **OFF時の後方互換**: `fluorescence: false` で従来と同一seed同一結果
     （既存の再現性テストのシーンで実装前後のハッシュ比較…は実装前の値を
     取れないので、「OFF時は `sample_fluorescence` が呼ばれない」ことの確認と
     既存テストが無傷であることで代替）

**完了条件**: 新規・既存テスト全通過。`chatcarlo run examples/chest_room.yaml
-n 1e5` が正常終了し、吸収/脱出割合が実装前とオーダーで変わらない
（chest_roomは150 kVp未満・鉛は薄いので蛍光の影響は小さいはず。
大きく変わったらバグを疑う）。

### Phase 3 — 数値検証と文書化

**タスク**:
1. **collision vs track-length の整合確認**: 鉛スラブを含むシーンで
   `--dose-grid` を実行し、鉛ボクセルの collision（energy_deposited系）と
   kerma track-length の比が蛍光ONで改善する（OFFより1に近づく）ことを
   確認、数値を記録。
2. **物理量の独立検算**: 100 keV→鉛の蛍光発生率の期待値を
   「K殻分率 × ω_K」から手計算し、シミュレーションの発生率と統計誤差内で
   一致することを確認（これは実装の自己整合でなく解析検算）。
3. `docs/lessons_learned.md` ではなく本計画書の末尾に結果サマリを追記し、
   CLAUDE.md の Physics 節に蛍光X線の1文（K殻のみ・カットオフ5 keV・
   `physics.fluorescence` トグル）を追記。README の機能一覧も同様。
4. コミット（コミットメッセージは日本語、既存スタイルに合わせる）。

**完了条件**: 検算一致・文書更新・コミット済み。

## 結果サマリ（実装完了、2026-07-17）

Phase 1〜3をすべて実施した。以下は実際に得られた数値。

**Phase 1データ層**: `fluorescence_k_data(82)`（鉛）で K端 88.00 keV・ω_K
0.9634・Kα1(KL3) 74.97 keV を確認（xraylib直接値と一致）。`fluorescence_k_data(74)`
（タングステン）で K端 69.53 keV・Kα1 59.32 keVを確認。`tests/test_fluorescence.py`
で全てスポットチェック済み。

**実装上の落とし穴**: `xraylib.CS_Photo_Partial(Z, K_SHELL, E)` は軽元素・
高エネルギー（例: Z=8酸素、E>~103 keV）で "Spline extrapolation is not
allowed" の `ValueError` を投げる。水を含むほぼ全シーンで即座に再現する
致命的な地雷だった。対策: 元素のK線が全て`_FLUOR_CUTOFF_KEV`(5 keV)未満
なら（軽元素はほぼ全てこれに該当— 酸素Kα~0.5 keV、カルシウムKα~3.7 keVも
未満）その元素は蛍光を絶対発生させないため、`CS_Photo_Partial`呼び出し自体を
スキップする早期continueを追加（`physics.py`の`sample_fluorescence`）。
物理的に正しい上に、問題のあるZ・E領域を実際に回避できる。

**エネルギー保存**: 100 keV→鉛スラブ(0.5cm)、n=5000、蛍光ON/OFF両方で
`|deposited + escaped - incident| / incident < 1e-9`（実測は machine epsilon
オーダー ~3e-16）を確認。

**蛍光発生率の解析検算**: 100 keV単色光子が鉛(Z=82)で光電吸収された場合、
K殻イオン化確率×ω_K = 0.796109 × 0.9634 = 0.766971（解析値）に対し、
シミュレーション実測 0.766505（n=2,000,000、5σ=0.001495以内で一致）。
線エネルギー分布もRadRateの分岐比どおり（Kα1 74.97 keV: 51.8%、
Kα2 72.80 keV: 30.8%、Kβ1 84.94 keV: 11.5%、Kβ3 84.45 keV: 6.0%）。

**collision vs track-length estimator整合**: 100 keV単色ビーム→鉛スラブ3cm
（ほぼ全数光電吸収）、n=500,000で比較（track-length/collision比、1に近いほど
整合）:

| 条件 | collision estimator [MeV/history] | track-length estimator [MeV/history] | 比 |
|---|---|---|---|
| 蛍光OFF | 0.09960 | 0.03765 | 0.378 |
| 蛍光ON | 0.08572 | 0.07863 | 0.917 |

蛍光実装により整合が大幅に改善（0.378→0.917）。残る差はコンプトン散乱後の
多重相互作用・境界効果等によるもので、単一の物理過程に起因するものではない
ため、この程度の残差は許容範囲と判断（完全な1.0はそもそも期待しない —
collision推定量とkerma推定量は定義上厳密に一致する量ではない）。

**scene.yaml統合**: `physics.fluorescence: true/false`をトグルとして追加、
デフォルトtrue。CLI (`chatcarlo run`) に「蛍光X線放出イベント数」を追加出力。
`chatcarlo trace`のtrajectory記録・preview.py凡例にも"fluorescence"イベント
（■マーカー）を追加し、正常に動作することを確認。

**テスト結果**: 新規11件含め既存全148件（実装完了時点）が通過
（`.venv/bin/python -m pytest tests/ -q`）。

**Phase 4（EGS5相互検証）は未着手**。次のセッションで`vive-crosscheck`により
実施する。

### Phase 4 — EGS5相互検証（別セッション、vive-crosscheck で実施）

本計画書のスコープ外だが後続として登録:
100 keV単色ビーム＋鉛スラブ透過のケースを、EGS5側 IEDGFL=1（K/L蛍光ON）
で `vive-crosscheck` により突き合わせる。許容基準は crosscheck 側で事前登録。
ChatCarloはK殻のみなので、EGS5のL蛍光分が系統差として出る可能性を
事前登録文書に明記しておくこと（鉛L線は透過しないため透過線量比較なら
差は小さいはず）。

## 実行者（Sonnet）への注意

- 触ってよいファイル: `chatcarlo/materials.py`, `chatcarlo/physics.py`,
  `chatcarlo/transport.py`, `chatcarlo/scene.py`, `chatcarlo/trajectory.py`
  （イベント種追加が必要な場合のみ）, `tests/test_fluorescence.py`（新規）,
  `CLAUDE.md`, `README.md`, 本計画書。**それ以外は変更しない。**
- tally.py・dose_coefficients.py は**触らない**。μen/ρは蛍光考慮済みの
  係数なので track-length 側に変更は不要（背景節参照）。
- 既存の設計パターン（`material_groups` でのグループ化、棄却法のpendingループ、
  `sample_*_element` の同型ヘルパ）を必ず踏襲する。新しい抽象は導入しない。
- 数値定数（K端エネルギー等）をコードやテストに手打ちで埋め込むのは
  テストの期待値（許容幅つき）のみ許可。実装本体は必ずxraylib経由。
- 各Phase完了時に `pytest tests/ -q` 全通過を確認してから次へ進む。
- 報告時は CLAUDE.md の報告ルールに従い、数値（蛍光発生率・エネルギー保存の
  残差・estimator比の改善）を本文に具体的に書くこと。
