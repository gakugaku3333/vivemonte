# 事前登録: 束縛コンプトン散乱（S(q)非干渉性散乱関数）実装

作成日: 2026-07-17。Phase 2b（[PDD_RESULTS.md](PDD_RESULTS.md)）で特定した残差の
うち、μen/ρ補間アーティファクト修正後も残る約−0.91%（対 標準EGS5 INCOH=0）の
主成分＝コンプトン散乱のエネルギー移行モデル差への対応。

## 背景・動機

現状のChatCarloは内部で不整合な状態にある: 相互作用種別の抽選
（`transport.py`の`p_compt`）は`xraylib.CS_Compt`（EPDLベース、**S(q)束縛補正込み**の
非干渉性散乱断面積）を使う一方、実際の散乱角・エネルギー抽選
（`physics.py`の`sample_klein_nishina`）は**自由電子Klein-Nishina**（束縛効果なし）。
つまり「コンプトンが起きる確率」は束縛込み、「起きたあとの角度・エネルギー移行」は
束縛なしという設計不整合がある。S(q)補正込みのサンプラーへの置き換えは、単なる
精度向上ではなくこの不整合の解消でもある。

## スコープ

- **対象**: S(q)非干渉性散乱関数によるコンプトン散乱角の補正サンプリングのみ
  （`xraylib.SF_Compt(Z, q)`、既存の`rayleigh_form_factor_table`と同じEPDLベース
  データソース）。
- **対象外（明示的に見送り）**: Doppler広がり（インパルス近似、束縛電子の運動量
  分布によるコンプトン端のぼけ）。理由: (1) 診断領域の線量計算への影響はS(q)角度
  補正よりさらに一桁小さいと予想される、(2) 比較対象のEGS5 INCOH=1もDoppler無し
  設定（IBOUND=1×INCOH=1はS(q)補正のみでDoppler広がりは含まない）であり、
  スコープを合わせないと相互検証の対象がずれる。
- **対象外**: 蛍光X線（特性X線）の放出・輸送。低Z材質では既存の判断（無視できる
  水準）を維持。

## 実装方針

Rayleigh散乱サンプリング（`sample_rayleigh_element`/`sample_rayleigh_cos_theta`/
`rayleigh_form_factor_table`）と対称的な設計にする:

1. `materials.py`: `compton_element_weights(material, energies_keV)` —
   `rayleigh_element_weights`と同型、`CS_Compt`版。化合物中でコンプトンがどの
   構成元素で起きたかを質量分率×元素別コンプトン断面積で抽選。
2. `materials.py`: `incoherent_sq_table(z)` — `rayleigh_form_factor_table`と同型、
   `xraylib.SF_Compt`版、lru_cacheでq∈[0,20]Å⁻¹をテーブル化。
3. `physics.py`: `sample_compton_bound(materials, energies_keV, rng)` —
   (a) 元素Zを`compton_element_weights`で抽選、(b) 既存の自由電子Kahn型棄却法
   （`sample_klein_nishina`のロジック）で(ε, cosθ)を提案、(c) S(Z,q)/Z
   （q=E·sin(θ/2)/hc）を追加の受理確率として棄却法を重ねる。
   S(Z,q)はq→∞でZに単調収束し常にS(Z,q)≤Zが成り立つため、
   S(Z,q)/Z∈[0,1]がそのまま有効な受理確率になる（追加の包絡線調整は不要）。
   既存`sample_klein_nishina`は削除しない（検証・比較用に残す）。
4. `transport.py`: コンプトン分岐の呼び出しを`sample_compton_bound`に切替。
   フラグは設けない単純切替（研究コードで新旧2系統を維持する保守コストを避ける）。
   旧自由電子KN挙動が必要な検証は`sample_klein_nishina`の直接呼び出しで可能。

## 判定基準（事前固定）

1. **サンプラー単体検証**: `sample_compton_bound`から得られる平均エネルギー
   移行率 <T>/E0（60 keV水、大標本）が、既存の独立机上検算
   `check_compton_transfer.py`の`average_transfer_fraction(60.0, use_sq=True)`
   が返す`t_Sq`と統計誤差内（標本標準誤差の3σ以内）で一致すること。
   一致しない場合は実装バグの可能性が高く、Phase 2b再検証に進まない。
2. **Phase 2b再検証（対EGS5 INCOH=1）**: `run_chatcarlo_pdd60.py`をS(q)版で
   再実行し、EGS5 INCOH=1再実行結果（既存、`pdd60_phantom_incoh1/egs5job.out`）
   との47ビン平均相対差が**±0.3%以内**に収まること。現状（PCHIP修正後、
   自由電子KN版）で+0.172%なので、これがさらに0に近づくか同水準を維持すると
   予測する。
3. **Phase 2b参考値（対EGS5 INCOH=0、予測を事前登録）**: 標準EGS5(INCOH=0)は
   角度分布は自由電子KNのままエネルギー付与のみ実測するため、S(q)化の影響は
   限定的と予想する。**予測: 47ビン平均相対差は現状の−0.91%から大きくは変化
   しない（±0.2pt程度の範囲に留まる）**。もし大きく変化した場合（例: 1pt以上
   動いた場合）は、S(q)化がEGS5 INCOH=0との比較で持つ意味を想定と異なる形で
   説明できていないことを意味し、原因究明を優先し「一致した」とは報告しない。
4. **既存回帰テスト**: `tests/test_transport.py`の一次透過率Beer-Lambert照合、
   `tests/test_tally.py`のcollision/track-length二推定量交差検証が、乱数消費列の
   変化を踏まえてもなお統計許容内で通ること（許容幅を超えて崩れた場合は
   実装バグを疑う）。

## 既知のリスク・留意点

- コンプトン分岐の乱数消費列がS(q)棄却法の追加により変わるため、**既存の
  全シミュレーション結果（`examples/`配下の期待値、他の相互検証成果物）は
  同一seedでも数値が変わる**。回帰テストで統計許容幅を使っている箇所は
  影響を受けないはずだが、厳密一致を要求しているテストがあれば洗い出す。
- 低エネルギー・重元素では S(Z,q)/Z の受理確率が小さくなる場面があり、
  棄却法の効率（採択率）が自由電子KN単体より低下する可能性がある。
  性能への影響は許容範囲か実装後に確認する（診断領域では極端に悪化しない
  はずだが、鉛遮蔽など重元素シーンで顕著な場合は記録する）。
- `SF_Compt`は`FF_Rayl`と同じEPDLデータソースだが、q=0近傍でxraylibが
  外挿エラーを返す既知の挙動（`check_compton_transfer.py`のクリップ処理と
  同様の対応が必要）。

## 完了後の文書化

- `PDD_RESULTS.md`: Phase 2b再検証結果を追記
- `lessons_learned.md`: 設計不整合（断面積は束縛込み・角度分布は自由電子）が
  存在していたこと自体を教訓として記録
- `CLAUDE.md`: 「物理」節のコンプトン記述をS(q)補正込みに更新
- vive-auditorによる監査（プロジェクト慣習）
