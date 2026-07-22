# 計画: 統計不確かさの一級市民化 — ボクセル相対誤差マップの標準搭載

作成日: 2026-07-22 / ステータス: **Phase 0・1・2 実装・検証完了。Phase 3（CLI/plot出力）
以降は未着手** / 実行担当: Claude

親コンテキスト: [[future-directions]]の優先候補2（「1の高速化 → 2の統計の見える化 →
3のタリー精密化」の2番目）。候補1（高速化）は`docs/plan_transport_speedup.md`・
`docs/plan_phase3_parallel.md`・`docs/plan_rayleigh_compton_importance_sampling.md`で完了済み。

この計画書は設計判断と実施順序を事前に固定してある。実行者は各ステップを上から順に
実施し、**受入基準を通してから次へ進む**こと。判断に迷う箇所は勝手に設計変更せず
「設計判断（確定事項）」に戻るか、ユーザーに確認する。

## 背景と目的

現状、統計誤差は**検証スクリプト側で毎回手書き**されている:

- `docs/egs5_crosscheck/run_chatcarlo_pdd60.py` は`per_history`配列を全history分
  確保して`arr.std(ddof=1)/np.sqrt(N)`を計算する。
- `docs/egs5_crosscheck/run_chatcarlo_water60.py` は二項近似
  `sqrt(p(1-p)/N)`を手計算する。
- `chatcarlo run --dose-grid`の出力には**統計誤差が一切出ない**。最大吸収線量・
  最大H\*(10)も、材料別吸収エネルギーも、点推定値だけが印字される。

これは3つの実害を生んでいる。

1. **「最大値」が信用できるか判定できない。** CLAUDE.mdの"Known sharp edges"に
   ある「解像度を細かくするほど最大値が増大する（max_substepsクランプ由来の
   極値統計と推定・未検証）」は、まさに統計誤差が見えないから未検証のままに
   なっている問題。最大値ボクセルの相対誤差Rが並記されれば、「Rが0.4の最大値」は
   その場で信用できないと分かる。
2. **遮蔽評価で「統計不足」を検知できない。** 監査で実際に出た所見
   （操作室遮蔽・鉛2cm背後でヒット数ゼロ、統計不足で評価不能。[[project-status]]）
   は、人間の監査官が気付いたから発覚した。相対誤差マップがあれば自動で分かる。
3. **分散低減（将来の候補）の効果を測る土台がない。** 分散低減の良し悪しは
   「同じ壁時間でRがどれだけ下がったか」でしか評価できない。Rを出す仕組みが
   先に要る。

**目的**: `run --dose-grid`にボクセルごとの相対誤差マップを標準搭載し、
「どこの統計が足りないか」を見えるようにする。あわせて`--dose-grid`なしでも
出る材料別吸収エネルギー等のスカラー量にも同じ推定器で誤差を付ける。

### 目標（事前登録）

- **主目標**: `run --dose-grid`が、最大吸収線量・最大H\*(10)に**相対誤差Rと
  寄与バッチ数を並記**し、`--dose-out`の.npzに**Rマップを保存**し、
  `chatcarlo plot`が**Rマップを描画**できる。
- **絶対制約**: 統計機構のON/OFFで**輸送結果（kerma_keV・h10配列を含む）が
  ビット一致**すること。乱数は一切消費しない。既存全テスト（233件）が通る。
- **副目標**: 相対誤差の推定器を`chatcarlo/tally.py`の公開純関数として提供し、
  EGS5相互検証スクリプト群が自前のΣx・Σx²計算を捨てられる状態にする。

## 設計判断（確定事項）

### 1. バッチ統計方式を採る（history-by-history方式は採らない）

MCNP流のhistory-by-history分散（各historyのボクセル寄与を個別に持つ）は、
本コードの構造では成立しない。輸送は`transport_photons`が**バッチ内の全光子を
同時にベクトル化して飛ばす**ため、`accumulate_track_length`のスコアリング時点で
「どのhistoryの寄与か」を分離して保持するには(photon_id, voxel)の疎ペアを
バッチ全体で溜める必要があり、代表条件（20万光子×平均数区間×最大40サブステップ）で
数千万要素・数百MBの一時領域になる。

採るのは**バッチ統計**: バッチ境界＝history境界なので、バッチ和 S_b は独立な
whole-history寄与の和である。したがってバッチ統計は**history単位の分散を
そのまま不偏推定する**（1光子が複数区間・複数ボクセルにまたがる相関は x_i の
内側に畳み込まれており、区間単位の分散推定が取り違える相関をここでは取り違えない）。
副次的に、タリーの層化乱数点によるサンプリング分散も同じRに含まれる——
つまりRは「報告値の総合的な統計不確かさ」であって物理サンプリング分だけではない。

### 2. 推定器: 不等バッチサイズでも不偏な形を使う

バッチ b の history数を n_b、そのバッチのボクセル寄与和を S_b、
総history数 N = Σn_b、バッチ数 M とする。

```
T = Σ S_b                      （既存のグリッド配列そのもの）
Q = Σ S_b² / n_b               （新規アキュムレータ）
σ̂² = (Q − T²/N) / (M − 1)      （1history あたり寄与 x_i の分散の不偏推定）
SEM = σ̂ / √N
R   = SEM / (T/N)  =  √(σ̂²) / (√N · T/N)
```

不偏性: E[Q] = Mσ² + Nμ²、E[T²/N] = σ² + Nμ²、差 = (M−1)σ²。**n_b が
バッチごとに違っても厳密に不偏**（証明は上記の期待値計算がn_bに依存しないこと
から従う）。この性質が効くのは、並列時にワーカーごとの端数バッチが必ず出るため。

**並列集約は (T, Q, M, N) の単純加算**で済む（既存のワーカー番号順の加算順序固定を
そのまま踏襲する）。

**分散低減への前方互換**: 将来ロシアンルーレット/スプリッティング等で光子重みが
1でなくなっても、S_b が重み付き寄与の和になるだけで推定器の形は変わらない
（historyでバッチを切っているため）。逆にイベント単位で数える方式は重み導入で
壊れる。同じ理由から、後述のヒット数マップも**イベント数ではなくバッチ数**で数える。

### 3. 統計バッチ ≡ 輸送バッチ。batch_size は自動調整しない

「M=30本になるようbatch_sizeを縮める」ことは**しない**。`batch_size`は
`rng.random(n)`の消費ブロック長を決めており、変えると同一seedでも光子ごとの
乱数が変わってビット一致が崩れる（`tests/test_reproducibility.py`と
`plan_phase3_parallel.md`が守っている不変条件）。

したがって **M = ceil(n_histories / batch_size)**（並列時は各ワーカーで同様、
総和がM）。既定batch_size=200,000では:

| n_histories | M | Rの相対不確かさ 1/√(2(M−1)) |
|---|---|---|
| 1e5 | 1 | 推定不能（Rはnan） |
| 1e6 | 5 | 35% |
| 1e7 | 50 | 10% |
| 1e8 | 500 | 3% |

**遮蔽評価スケール（n≥1e7、この機能が本当に要る領域）では十分**。小nでM<2の
ときは黙って嘘の数字を出さず「バッチ数不足で誤差推定不可」と明示する。
より細かい制御が要るユーザーのために`run --batch-size`をCLIに露出するが、
**既定値は変えない**（変えると結果がビット一致しなくなる旨をhelpに書く）。

**帰結（意識的に受け入れる）: 既定の`run --dose-grid`（`-n`既定1e5）ではRが出ない。**
M=1なので「誤差推定不可」メッセージだけが出る。これを避けるために`--dose-grid`時の
既定`-n`を引き上げることは**しない**——既定コマンドは「まず動かして絵を見る」
教育・下見用途で秒オーダーに収まっていることに価値があり、統計のために分オーダーへ
落とすのは割に合わない。代わりに**メッセージを実行可能な指示にする**:

```
  統計誤差: バッチ数が不足しています（M=1、必要なのはM>=2、実用上はM>=20）。
  現在の設定: n_histories=1e5 / batch_size=200,000 → M=ceil(1e5/200,000)=1
  対処: -n を 4e6 以上（M>=20）にするか、--batch-size を 5000 に下げてください。
```

数値は実際の設定から計算して埋める（固定文言にしない）。

**同一nでもserialとparallelでMが違う点**も注記する: n=1e6・batch_size=2e5なら
serialはM=5、workers=4では各ワーカーが自分の25万historyを切り上げて2バッチ持つため
M=8になる。推定器は不等バッチで不偏なので値は正しいが、**Rの数値は目に見えて変わる**
（`--workers`を変えるとビット一致しない・統計的同等のみ、という既存の契約の延長）。
バグ報告と誤解されないようhelpとREADMEに1行書く。

### 4. 実装は「スナップショット差分」方式（スクラッチ配列方式は採らない）

バッチ寄与 S_b の取り出し方は2案ある:

- (a) スクラッチ配列にスコアし、バッチ末に total へ畳み込む
- (b) **これまで通り total へ直接スコアし、バッチ末に delta = total − prev、
  prev ← total とする（採用）**

(b)を採る理由:
- **totalの浮動小数点加算順序が現行と完全に同一**なので、統計機構ON/OFFで
  kerma_keV/h10配列が**ビット一致**する（絶対制約）。(a)は畳み込みの分だけ
  加算順序が変わり最終ビットがずれる。
- **失敗モードが安全側**。`transport_photons`をテスト・検証スクリプトから直接
  呼ぶ既存利用者がバッチ終端処理を呼ばなくても、(b)なら総和は正しいまま
  （統計量だけが出ないだけ）。(a)だと総和がゼロのままという致命的な沈黙バグになる。
- 差分の桁落ちは無視できる: delta の絶対誤差は ~eps·T ≈ M·eps·S_b、相対で
  約1e-14。真に寄与ゼロのバッチでは total がビット単位で不変なので delta は
  厳密に0になる（偽の微小値は出ない）。

### 5. 相対誤差はカーマ/H\*(10)の生アキュムレータから計算する

吸収線量[Gy] = カーマ/質量、H\*(10)[pSv] = 飛程積分/体積 は、いずれもボクセルごとに
**決定的な定数倍**である。したがって R_dose ≡ R_kerma、R_h10 ≡ R_h10track。
校正係数（mAs基準の実光子数）も決定的な定数倍なので、**Rマップは相対値出力にも
絶対値出力にも同じものが使える**。

例外の注意書きを1行入れること: `ctdi_vol_mGy`校正だけは校正係数自体が別のMC計算
（`ctdi_per_history_Gy`）由来なので、絶対値の総合不確かさはRだけでは足りない。
その定量化は本計画の非目標（下記）。

### 6. 既定でON、逃げ道は `--no-uncertainty`

「標準搭載」なので`--dose-grid`指定時は既定でONにする。代償はメモリで、
ボクセルあたり 16バイト → 約52バイト（kerma: total/prev/Q ＝ 8×3、h10: 同、
共通ヒット数 int32 ＝ 4）の**約3.25倍**。`__main__.py`の並列メモリ警告の
概算式（現行 `prod(shape)*8*2*(workers+1)`）を更新する。細解像度×多ワーカーで
足りない場合の逃げ道として`--no-uncertainty`を用意する。

ヒット数マップを kerma/h10 で共有してよい理由: μen/ρ も h\*(10)/Φ も
全エネルギー・全材料で正であるため、ある区間が触れるボクセル集合は両者で同一。

### 7. Rの解釈ガイドを併せて出す（Rは必要条件であって十分条件ではない）

MCNP流の目安（R<0.05: 一般に信頼できる / 0.05–0.10: おおむね信頼できる /
0.10–0.20: 疑わしい / >0.20: 意味を持たない）をREADMEとCLI出力の注記に載せる。

**同時に、この機能の目玉である「鉛2cm背後」ケースこそRが嘘をつくことを明記する**:
数本のhistoryしか届いていないボクセルでは寄与分布がゼロ膨張・強い歪みを持ち、
「まだ大きな寄与を引いていないだけ」の状態でRが小さめに出て**偽の安心**を与える。
これはMCNPでも古典的な既知の落とし穴。だから**Rマップ単独では出さず、寄与バッチ数
マップを必ず併記し、寄与バッチ数が少ないボクセルのRは信用しない**という運用を
ドキュメント・CLI警告の両方で強制する。

## 非目標（やらないこと）

- **分散低減そのもの**（スプリッティング/ロシアンルーレット/強制衝突）。本計画は
  その効果を測る土台を作るだけ。
- **history-by-history分散**（設計判断1）。
- **CTDI校正係数自体の不確かさ伝播**（設計判断5の例外）。注記のみ。
- **FOM（= 1/(R²·T)）の実装**。基準ボクセルの定義を決める必要があり、分散低減に
  着手するときに一緒に決めるほうが手戻りがない。申し送りに書く。
- **`accumulate_track_length`のサブステップ方式そのものの変更**（[[future-directions]]
  候補3「解析的重なり長方式」は別計画）。ただし本計画の成果は候補3の検証手段に
  なる（解像度を細かくしたとき最大値の増大が「Rの範囲内の揺らぎ」か「系統的増大」か
  を初めて切り分けられる）。

## フェーズ計画

### Phase 0 — 事前計測（コード変更なし）

1. 代表2条件でベースラインの壁時間とピークメモリを測る:
   - `examples/chest_room.yaml`、n=1e6、`--dose-grid --resolution 5`、workers=1
   - 同、n=1e7、workers=4
2. `docs/speedup_baseline/`に既存の計測ファイル群と同じ形式で記録する。
3. **受入**: 数値が記録されている（Phase 4のオーバーヘッド判定の基準になる）。

### Phase 1 — 推定器とアキュムレータ（tally.py、輸送非依存）

1. `VoxelGrid`に`track_uncertainty: bool = False`（既定OFF＝既存利用者に非侵襲）と、
   ONのとき確保する配列 `kerma_sum2` / `h10_sum2` / `_kerma_prev` / `_h10_prev` /
   `n_batches_hit`(int32) を追加。`n_batches` / `n_histories` のスカラーも持つ。
2. `VoxelGrid.end_batch(n_histories_in_batch: int)` を追加（設計判断4のスナップショット
   差分）。`track_uncertainty=False`なら何もしない。
3. 純関数として公開する:
   - `relative_error(total, sum2, n_batches, n_histories) -> np.ndarray`
     （total==0 のボクセルは nan、M<2 は全nan、σ̂²の丸め負値は0にクランプ）
   - `standard_error(...)`（同上のSEM。スカラーにもグリッドにも使える形にする）
   - `combine_moments(...)`（並列集約の (T,Q,M,N) 加算。純関数にしてテストしやすく）
4. スカラー版アキュムレータ`ScalarMoments`（材料別吸収エネルギー・吸収/脱出割合用。
   dict of 材料名 -> (S, Q) を持ち、あるバッチに現れない材料は S_b=0 として扱う）。
5. **テスト（`tests/test_uncertainty.py`新規）**:
   - **代数の厳密検証（batch_size=1）**: 既知の乱数列から作った合成history列を
     1historyごとに1バッチとして統計を取ると、Q=ΣS_b²/n_bはn_b=1のとき
     Σx_i²に厳密に一致するため、σ̂²は`x.std(ddof=1)**2`と数式的に完全一致する。
     `relative_error`が`x.std(ddof=1)/sqrt(N)/mean`と`rtol=1e-12`で一致することを
     確認する。
   - **不等バッチサイズでの不偏性（統計的検証、厳密一致ではない）**: バッチサイズを
     不揃いにすると、S_b²/n_bはΣ_{i∈b}x_i²の代わりにはならない（交差項が残る）ため
     σ̂²はbatch_size=1の値と数式的には一致しない——**別の不偏推定量になるだけ**
     （設計判断2の期待値計算はn_bに依存せず成り立つ）。既知分布からの乱数を
     不等バッチに区切り、多数回反復した σ̂² の平均が真の分散に統計的に収束する
     ことを確認する（許容誤差付きの統計検証であり、単発の厳密一致は要求しない）。
   - `combine_moments`: バッチ列を2分割して個別集計→合成した結果が、
     一括集計と`rtol=1e-12`で一致（こちらは単純加算の結合律なので厳密一致）。
   - 退化ケース: M=1でnan・例外なし、mean=0でnan（infでも警告でもない）、
     σ̂²が丸めで微小負になる入力で0にクランプ。
6. **受入**: 上記テストが通り、既存テストに影響なし（この時点で輸送は未変更）。

### Phase 2 — 輸送への組み込み（transport.py）

1. `run_transport(..., track_uncertainty: bool = True)` を追加。`dose_grid=True`の
   ときだけ意味を持つ。`VoxelGrid`生成時に渡す。
2. `_run_batches`のループ末尾で`grid.end_batch(n)`と`ScalarMoments.end_batch(n)`を呼ぶ。
   **`transport_photons`は一切変更しない**（乱数消費・スコアリング経路を触らない）。
3. 並列パス: ワーカー戻り値に `kerma_sum2` / `h10_sum2` / `n_batches_hit` /
   `n_batches` / スカラー版モーメントを追加し、既存のワーカー番号順ループで加算する
   （加算順序固定の理由は`plan_phase3_parallel.md`設計判断のまま）。
4. `TransportResult`に統計フィールドを追加（`n_batches`、スカラー量のSEM、
   グリッドのRを返すアクセサ）。
5. **テスト**:
   - **ビット一致（最重要）**: 同一seed・同一workersで`track_uncertainty`
     True/Falseの2回runし、`grid.kerma_keV`と`grid.h10_track_pSv_cm3`が
     `np.array_equal`で完全一致。材料別吸収エネルギーも完全一致。
   - **1/√Nスケーリング**: n=4Nの run のRが n=N の run のR のおよそ1/2
     （代表ボクセルで20%以内）。
   - **バッチ分割不変性（統計的）**: 同一条件でbatch_sizeだけ変えた2runのRが
     数十%以内で一致（Rは点推定なので厳密一致は要求しない）。
   - **並列一致**: workers=1とworkers=4のRが同オーダー
     （`tests/test_parallel_transport.py`の既存スタイルに合わせる）。
   - **ブルートフォース照合**: 小さなシーン・小さなグリッドで、`batch_size=1`
     （＝S_bが1historyの寄与そのもの）で走らせたときのRが、同じ寄与列から
     直接計算した`std(ddof=1)/√N/mean`と`rtol=1e-10`で一致。
     **監査官が「二項近似ではなく直接std」と修正した教訓（[[future-directions]]）を
     ここで恒久化する**。
6. **受入**: 上記テスト＋既存233テストが全通過。

### Phase 3 — 出力層（__main__.py / diagnostics.py / plotting.py）

1. **`run`の出力**（この計画の目玉）:
   ```
     最大吸収線量 [Gy/history]: 1.234e-09  (相対誤差 R=0.043, 寄与バッチ 48/50)
     最大H*(10)   [pSv/history]: 5.678e-04  (相対誤差 R=0.112, 寄与バッチ 12/50)
     グリッド統計: R<0.10 のボクセル 63.2% / 寄与バッチ0（未到達）のボクセル 11.4%
   ```
   - 材料別吸収エネルギーにも `± SEM (相対x.xx%)` を付ける（`--dose-grid`不要）。
   - M<2 のときは値の代わりに「バッチ数不足（n_histories/batch_size < 2）で
     誤差推定不可」と出す。
2. **診断警告（diagnostics.py）**: 既存の`background_medium_warning` /
   `near_source_air_warning`と同じ形で純関数を足す:
   - `unreliable_max_warning(R, n_hit_batches, n_batches)` — 最大値ボクセルの
     R>0.10、または寄与バッチ数が全体の1/4未満のとき警告。文面には設計判断7の
     「寄与バッチ数が少ないとRは楽観側に嘘をつく」を含める。
   - 既存の"Known sharp edges"（解像度を細かくすると最大値が増大する件）に、
     Rを見て判断せよという誘導を入れる。
3. **npz出力**: `rel_err_dose` / `rel_err_h10` / `sem_dose_per_history_Gy` /
   `sem_h10_per_history_pSv` / `n_batches` / `n_batches_hit` を追加保存
   （`--no-uncertainty`時は書かない）。既存キーは変更しない（後方互換）。
4. **`plot`**: `--quantity` に `relerr-dose` / `relerr-h10` を追加。
   - LogNormではなく線形（0〜0.5程度でクリップ）、`inferno`とは別のカラーマップ、
     mean=0のボクセルはマスク。
   - 寄与バッチ数が閾値未満のボクセルは別色でマスクし、凡例で「統計未到達」と示す
     （設計判断7の運用強制）。
   - 旧npz（新キーなし）を渡されたら親切なエラーで落とす。
5. **受入**: CLIスモークテスト（npzに新キーが入る／`--no-uncertainty`で入らない）、
   `plot --quantity relerr-dose`のスモークテスト、`unreliable_max_warning`の単体テスト。

### Phase 4 — 実測・検証・文書化

1. **オーバーヘッド実測**（Phase 0のベースラインと比較）:
   - 壁時間の増分が**+10%以内**なら追加最適化しない。
   - 超えた場合のみ、`accumulate_track_length`が既に持っている`flat_idx`から
     「そのバッチで触れたボクセル」だけを`np.unique`で拾い、end_batchの
     O(M×n_vox)全面走査を触れたボクセルだけに絞る最適化を検討する。
     **推測で先に実装しない**（`docs/lessons_learned.md`「高速化の対象は推測せず
     実測してから決める」）。
   - メモリ実測値で`__main__.py`の警告概算式（設計判断6）を裏取りする。
2. **動く実例を作る**（CLAUDE.mdの報告ルール: 成果物を提示する）:
   - `examples/chest_room.yaml`のRマップPNG。
   - **遮蔽背後ケース**: 鉛遮蔽を含むシーンで、遮蔽背後がRマップ上で
     「寄与バッチ0」として塗り分けられる図。この機能の存在理由そのもの。
   - n=1e6 と n=1e7 のR比較（1/√Nで下がることを図で見せる）。
3. **検証スクリプトの置き換え**: `docs/egs5_crosscheck/run_chatcarlo_pdd60.py`の
   手書きSEM計算を`chatcarlo.tally`の公開関数に置き換え、**RESULTS.mdに記載済みの
   数値が変わらない**ことを確認する（副目標の実証。数値が動いたら推定器のバグ）。
4. **文書更新**:
   - `README.md`: Rの解釈ガイド（設計判断7）と「Rは必要条件であって十分条件ではない」。
   - `CLAUDE.md`: "Known sharp edges"にRによる判定方法を追記、コマンド例を更新。
   - `docs/lessons_learned.md`: 得られた教訓。
   - この計画書に実施記録を追記。
   - メモリ更新: [[future-directions]]の候補2を完了に、[[project-status]]に現況。

## 受入基準（まとめ）

1. 統計ON/OFFで`kerma_keV`・`h10_track_pSv_cm3`・材料別吸収エネルギーが**ビット一致**。
2. `batch_size=1`でのRが直接計算の`std(ddof=1)/√N/mean`と`rtol=1e-10`で一致。
3. 不等バッチサイズでの`combine_moments`が一括集計と`rtol=1e-12`で一致。
4. Rが1/√Nで下がる（4倍history で20%以内の精度で半減）。
5. 既存233テスト＋新規テストが全通過。
6. `run --dose-grid`が最大値にRと寄与バッチ数を並記し、統計不足時に警告を出す。
7. `plot --quantity relerr-dose`が図を出し、統計未到達ボクセルを塗り分ける。
8. 壁時間オーバーヘッド+10%以内（超過時はPhase 4-1の最適化を実施してから再判定）。
9. `run_chatcarlo_pdd60.py`を公開関数に置き換えてもRESULTS.mdの数値が変わらない。

## リスクと落とし穴

- **M=1で無言のnan**。既定n_histories=1e5・batch_size=200,000ではM=1になる
  （設計判断3の「帰結」）。nanを黙って印字せず、必ず実際の設定値と対処法を
  埋め込んだメッセージにする（Phase 3-1）。
- **「Rが小さい＝正しい」と読まれる**。設計判断7の通り、寄与バッチ数マップを
  必ずセットで出す。ドキュメントで「必要条件であって十分条件ではない」を明言。
- **メモリ3.25倍**。細解像度×多ワーカーで既存の4GB警告に容易に触れる。
  警告式の更新（設計判断6）を忘れない。
- **`transport_photons`を直接呼ぶ既存コードとの相互作用**。`track_uncertainty`
  既定OFF＋スナップショット差分方式（設計判断4）で総和は常に正しいが、
  「end_batchを呼ばないと統計だけ出ない」ことを`VoxelGrid`のdocstringに明記する。
- **並列の加算順序**。(T,Q,M,N)いずれもワーカー番号順の加算を崩さない
  （崩すと同一(seed,workers)の再現性が失われる）。
- **`plot`の後方互換**。旧npzを読んだときに`KeyError`で不親切に落ちないようにする。

## 申し送り（本計画の外）

- **FOM = 1/(R²·T)**: 分散低減に着手するときに、基準ボクセルの定義と併せて決める。
- **候補3（解析的重なり長タリー）との接続**: 本計画のRがあれば、
  「解像度を細かくすると最大値が増大する」現象が統計揺らぎか系統誤差かを
  初めて切り分けられる。候補3の着手時に最初にやる測定として推奨。
- **CTDI校正係数の不確かさ伝播**（設計判断5の例外）。

## 実施記録（2026-07-22）

### Phase 0: ベースライン計測

`docs/speedup_baseline/uncertainty_phase0_baseline.txt`に記録。
- 条件1（chest_room.yaml, n=1e6, --dose-grid --resolution 5, workers=1）:
  real 27.57s、ピークメモリ約1.80GiB、グリッド形状(100,100,79)。
- 条件2（同, n=1e7, workers=4）: real 82.34s。
- Phase 4のオーバーヘッド判定基準（+10%以内）として、それぞれ30.3s・90.6sを
  上限値に設定済み。

### Phase 1: 推定器とアキュムレータ（tally.py）

設計判断どおり実装:
- `VoxelGrid`に`track_uncertainty: bool = False`と`kerma_sum2` /
  `h10_sum2` / `n_batches_hit` / `n_batches` / `n_histories` /
  `_kerma_prev` / `_h10_prev`を追加。`from_bbox`に`track_uncertainty`引数を追加
  （既定False、既存呼び出し元は無変更で動く）。
- `VoxelGrid.end_batch(n_histories_in_batch)`をスナップショット差分方式で実装
  （設計判断4）。`track_uncertainty=False`のときは何もしない。
- 純関数 `standard_error` / `relative_error` / `combine_moments` を追加。
  スカラー・ndarrayの両方に同じ式で動く（`np.asarray`で内部吸収）。
  M<2・mean=0はnan、丸め負分散は0クランプ。
- `ScalarMoments`（材料別吸収エネルギー等の辞書量向け）を追加。
  未出現材料をS_b=0として扱う仕様を実装・テストで確認。
- `VoxelGrid.kerma_relative_error()` / `h10_relative_error()` を追加
  （Phase 3の出力層から呼ぶための薄いラッパー）。

### 検証結果

`tests/test_uncertainty.py`（新規16件）:
- batch_size=1でのブルートフォース一致（`rtol=1e-12`）— relative_error・
  standard_error双方。
- 不等バッチサイズでの不偏性を4000回反復のモンテカルロで統計的に確認
  （真の分散との相対誤差5%未満）。batch_size=1の値とは数式的に一致しない
  ことを設計判断どおりdocstring・テストコメントで明記（当初計画の文言の
  誤り——「不等バッチサイズでも一致する」を「不等バッチサイズでも不偏」に
  訂正済み、本文参照）。
- `combine_moments`のスカラー・ndarray両形式が一括集計と`rtol=1e-12`で一致。
- 退化ケース（M<2、mean=0、丸め負分散）がすべてnan/0クランプで例外なし。
- `VoxelGrid`: 統計ON/OFFで`kerma_keV`・`h10_track_pSv_cm3`が
  **`np.array_equal`で完全ビット一致**（絶対制約の直接確認）。
  `track_uncertainty=False`時の`end_batch`呼び出しがno-opであることも確認。

既存テストスイート: 262件（新規16件込み）全通過、リグレッションなし
（Phase 1は`transport.py`を一切変更していないため非侵襲性は設計上自明だが、
念のためフルスイートを実行して確認）。

### Phase 2: 輸送への組み込み（transport.py）

設計判断どおり実装、`transport_photons`は無変更のまま:
- `_run_batches`のループ末尾で毎バッチ`grid.end_batch(n)`（gridがNoneでなければ
  常に呼ぶ——grid.track_uncertainty=Falseならno-opなので条件分岐不要）と
  `energy_moments.add_batch(result.energy_deposited, n)`（track_uncertainty=True
  のときのみ）を呼ぶよう変更。戻り値dictに`kerma_sum2`/`h10_sum2`/
  `n_batches_hit`/`grid_n_batches`/`grid_n_histories`/`energy_moments`
  （`ScalarMoments.as_moments()`の出力）/`scalar_n_batches`/`scalar_n_histories`
  を追加（並列ワーカーからのpickle転送用、serial時は既存のgrid直接参照と
  実質重複するが構造を統一するため常に含める）。
- `ScalarMoments.merge_from(other_moments, other_n_batches, other_n_histories)`
  を新規追加（tally.py）。他プロセスの`as_moments()`出力を材料ごとに単純加算し、
  自分に無い材料キーは(0,0)として扱う（`add_batch`の並列版）。
- `_run_worker`に`track_uncertainty`引数を追加し、ワーカー自身の
  VoxelGrid/ScalarMomentsへ独立に積む。
- `run_transport`に`track_uncertainty: bool = True`（既定ON、設計判断6）を追加。
  serial/parallel問わず`energy_moments`をrun_transport側で組み立て直す
  （serialは`agg`の`energy_moments`をmerge_from、parallelはワーカー番号順の
  ループ内で逐次merge_from——既存の加算順序固定パターンをそのまま踏襲）。
  `TransportResult`に`n_batches`/`energy_deposited_sem_MeV`/
  `energy_deposited_rel_err`を追加。ボクセルRは`result.grid.kerma_relative_error()`/
  `h10_relative_error()`を直接呼ぶ形（重複保持しない）。

**帰結**: `dose_grid=False`でも`track_uncertainty=True`（既定）なら材料別
付与エネルギーのSEM/Rが常に得られる（副目標の実現、grid不要）。

### 検証結果（Phase 2）

`tests/test_uncertainty_transport.py`（新規11件）:
- **ビット一致（最重要）**: serial・parallelともtrack_uncertainty ON/OFFで
  `kerma_keV`・`h10_track_pSv_cm3`・`energy_deposited_MeV`・
  `fraction_absorbed`/`fraction_escaped`/`mean_scatter_events`・
  `n_fluorescence`・`n_photons_real`が完全一致。
- **並列再現性の統計版**: 同一(seed,workers=2)の2回実行で`kerma_sum2`・
  `h10_sum2`・`n_batches_hit`・`n_batches`・材料別R まで完全一致
  （既存のplan_phase3_parallel.mdの契約の統計量への拡張）。
- **1/√Nスケーリング**: 単発runでは統計揺らぎが大きい（実測でratio
  0.98〜3.4と広く散ることを事前に確認）ため、6回反復平均で検証。
  4倍historyでRがおよそ半分（許容 1.3〜3.0倍、実測1.87倍）。
- **バッチ分割の統計的不変性**: 同一n_historiesでbatch_size違い
  （M=20 vs M=10）でもRは同オーダー（許容 0.2〜5.0倍、実測比1.45倍）。
- **並列/直列の同オーダー**: workers=1と2でRは同オーダー
  （実測比1.4倍程度、許容0.2〜5.0倍）。
- **ブルートフォース照合（最重要の一つ）**: `run_transport`を介さず低レベル
  API（`transport_photons`+`VoxelGrid.end_batch`を1historyずつ直接呼ぶ）で
  400historyを実際の物理輸送し、記録したper-history寄与列の直接計算
  `x.std(ddof=1)/√N/mean`と`grid.kerma_relative_error()`が`rtol=1e-9`で一致
  ——Phase 1で検証済みの推定器代数が、実際の輸送物理を通しても正しく
  配線されていることを確認。
- 副目標: `dose_grid=False`でも材料別SEM/Rが得られること、
  `track_uncertainty=False`でこれらが空dictになること、
  M=1（バッチ数不足）でR全体がnanになることを確認。

既存テストスイート: 273件（Phase 1・2の新規27件込み）全通過、
リグレッションなし。CLIスモークテスト（`chatcarlo run --dose-grid`、
既定track_uncertainty=Trueが暗黙に効く状態）も出力に変化なく正常動作を確認
（Phase 3でR/SEMを出力に追加するまでは既存出力のまま）。

### Phase 2完了後の批判的セルフレビュー（同日、ミューテーション試験）

「テストが通った」ではなく「**壊したらテストが落ちるか**」を確認するため、
推定器の各項を意図的に壊してテストを走らせた（4パターン）。

| ミューテーション | 検出 |
|---|---|
| `ScalarMoments.add_batch`の`/n_b`を削除 | ✅ 1件失敗 |
| `_batch_variance`の`M-1`→`M` | ✅ 4件失敗 |
| `standard_error`の`/√N`を削除 | ✅ 5件失敗 |
| **`VoxelGrid.end_batch`の`/n_b`を削除** | ❌ **27件全通過（検出できず）** |

**原因**: グリッドのRの絶対値を検証していたテストが全て n_b=1
（`/n_b`が恒等変換になる条件）だった——Phase 2のブルートフォース照合テストは
1historyずつ`end_batch(1)`を呼ぶ設計、Phase 1のグリッドテストは有限性しか
見ていなかった。**本計画で最も新規性が高く、不等バッチを成立させている当の
正規化が無検証だった**（`ScalarMoments`側は n_b=1000 と n_b=1 を突き合わせる
テストを偶然書いていたため検出された）。

**対処（テスト2件追加）**:
- `test_end_batch_normalizes_by_batch_history_count`: n_b をバッチごとに
  変えて（4,1,7,2,5,3）決め打ち値を積み、実装の`Q=ΣS_b²/n_b`とは**独立な式**
  ——分散分析の群間平方和 SSB=Σn_b(x̄_b−x̄)²——から期待値を組んで`rtol=1e-12`で照合。
- `test_relative_error_invariant_to_batch_grouping`: 同一のper-history寄与列を
  n_b=1 と n_b=10 でグループ化し、σ̂²が同じ母分散を推定していることを確認
  （正規化が抜けていれば約√10倍ずれる）。

両テストとも、修正後に同じミューテーションを再投入して**確かに落ちること**を
確認済み（テストがテストとして機能していることの確認）。

### 同レビューで発見・修正した実装上の問題

1. **メモリ警告式の未更新（実害あり）**: `run_transport`の既定を
   `track_uncertainty=True`にした一方、`__main__.py`のメモリ概算式は
   `8*2*(workers+1)`（＝統計なしの16バイト/ボクセル）のままだった。
   chest_room 1cm解像度・workers=4で**実態25.7GBを7.9GBと報告**する、
   安全側でない誤りだった。52バイト/ボクセルへ更新し、`--no-uncertainty`の
   案内も警告文に含めた。
2. **逃げ道の欠落（計画の実施順序ミス）**: 設計判断6は「既定ON＋
   `--no-uncertainty`で逃げる」だったが、既定ONをPhase 2で入れ、フラグを
   Phase 3送りにしていたため、ユーザーが3.25倍のメモリを回避する手段が
   ない期間ができていた。`--no-uncertainty`をPhase 2へ前倒しして解消。
3. **親プロセスでのスナップショット死蔵**: `_kerma_prev`/`_h10_prev`を
   `__post_init__`で確保していたが、並列時の親グリッドは集約先であり
   `end_batch`を呼ばない。1cm解像度で16バイト/ボクセル＝約1.6GBが完全な
   死蔵になっていた。**最初の`end_batch`での遅延確保**に変更。
4. **バッチごとの配列再確保**: `self._kerma_prev = self.kerma_keV.copy()` は
   毎バッチ グリッド2枚分を確保・解放し直す（細解像度では1枚790MB級）。
   `np.copyto`で既存バッファ再利用に変更。

修正後: 全275件通過。CLIで`--no-uncertainty`有無の最大吸収線量が一致
（ビット一致の実地確認）。実測オーバーヘッドは chest_room n=1e6
`--dose-grid --resolution 5` workers=1 で 27.57s→27.67s（**+0.4%**、
Phase 4の受入基準+10%に対し十分小さい）。

### 次のステップ

Phase 3（CLI出力: `run --dose-grid`の最大値へのR・寄与バッチ数の並記、
`diagnostics.py`への`unreliable_max_warning`追加、npz出力へのRマップ追加、
`plot --quantity relerr-dose/relerr-h10`の追加）が次の着手点。
なお`--no-uncertainty`フラグ自体は上記レビューでPhase 2へ前倒し済み。

**Phase 4へ申し送り（レビューで新たに認識した測定項目）**: 並列時に
ワーカーが親へ返す配列が2枚→5枚（16→36バイト/ボクセル）に増えたため、
pickleの転送量とピークメモリのスパイクが増える。Phase 4のメモリ実測では
壁時間だけでなく**ワーカーごとのRSSと集約時のスパイク**も採ること。
