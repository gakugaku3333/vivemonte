"""線源スペクトル生成（SpekPy統合）のテスト。"""
from __future__ import annotations

import numpy as np

from chatcarlo import spectrum


def test_spekpy_available():
    assert spectrum._HAS_SPEKPY, "spekpyがインストールされていません（.venv/bin/pip install spekpy）"


def test_spectrum_normalized_and_bounded():
    e, w = spectrum._default_spectrum(kvp=120.0, filtration_mm_al=2.5, anode_angle_deg=12.0)
    assert np.isclose(w.sum(), 1.0)
    assert np.all(w >= 0)
    assert e.min() > 0
    assert e.max() <= 120.0


def test_spectrum_mean_energy_is_plausible():
    """典型的な濾過条件では実効エネルギーはkVpの3〜6割程度になる（教科書的な目安）。"""
    e, w = spectrum._default_spectrum(kvp=120.0, filtration_mm_al=2.5, anode_angle_deg=12.0)
    mean_e = np.sum(e * w)
    assert 0.3 * 120.0 < mean_e < 0.6 * 120.0


def test_higher_filtration_hardens_spectrum():
    """濾過を厚くすると低エネルギー光子が優先的に除去され、平均エネルギーが上がるはず。"""
    e_thin, w_thin = spectrum._default_spectrum(kvp=100.0, filtration_mm_al=1.0, anode_angle_deg=12.0)
    e_thick, w_thick = spectrum._default_spectrum(kvp=100.0, filtration_mm_al=5.0, anode_angle_deg=12.0)
    mean_thin = np.sum(e_thin * w_thin)
    mean_thick = np.sum(e_thick * w_thick)
    assert mean_thick > mean_thin


def test_sample_spectrum_draws_within_range(monkeypatch):
    rng = np.random.default_rng(0)
    src = {"kvp": 80.0, "filtration_mm_al": 2.5, "anode_angle_deg": 12.0}
    energies = spectrum.sample_spectrum(src, 5000, rng)
    assert energies.min() > 0
    assert energies.max() <= 80.0


def test_explicit_spectrum_overrides_spekpy():
    rng = np.random.default_rng(0)
    src = {"kvp": 80.0, "spectrum": [{"energy_keV": 50.0, "weight": 1.0}]}
    energies = spectrum.sample_spectrum(src, 100, rng)
    assert np.all(energies == 50.0)


def test_fallback_spectrum_used_when_spekpy_unavailable(monkeypatch):
    monkeypatch.setattr(spectrum, "_HAS_SPEKPY", False)
    e, w = spectrum._default_spectrum(kvp=100.0, filtration_mm_al=2.5, anode_angle_deg=12.0)
    assert np.isclose(w.sum(), 1.0)
    assert e.max() <= 100.0


def test_export_import_caches_avoids_recomputation(monkeypatch):
    """export_caches/import_caches（並列ワーカーのSpekPy固定費削減用、
    docs/plan_phase3_parallel.md「積み残し」節）— import後は同じキーで
    SpekPyを再呼び出ししないこと。"""
    spectrum._spectrum_cache.clear()
    spectrum._heel_cache.clear()
    e1, w1 = spectrum._default_spectrum(kvp=90.0, filtration_mm_al=2.5, anode_angle_deg=12.0)
    spec_cache, heel_cache = spectrum.export_caches()
    assert (90.0, 2.5, 12.0) in spec_cache

    spectrum._spectrum_cache.clear()
    spectrum._heel_cache.clear()
    spectrum.import_caches(spec_cache, heel_cache)

    calls = []
    real_spek = spectrum._spekpy.Spek

    def _tracking_spek(*args, **kwargs):
        calls.append((args, kwargs))
        return real_spek(*args, **kwargs)
    monkeypatch.setattr(spectrum._spekpy, "Spek", _tracking_spek)

    e2, w2 = spectrum._default_spectrum(kvp=90.0, filtration_mm_al=2.5, anode_angle_deg=12.0)
    assert calls == []  # キャッシュヒットでSpekPyは一切呼ばれない
    assert np.array_equal(e1, e2) and np.array_equal(w1, w2)


def test_export_import_caches_is_a_pure_snapshot():
    """export_caches()は呼び出し時点のスナップショット（後続のキャッシュ更新は
    エクスポート済み辞書に反映されない）であること。"""
    spectrum._spectrum_cache.clear()
    spectrum._default_spectrum(kvp=70.0, filtration_mm_al=2.5, anode_angle_deg=12.0)
    spec_cache, _ = spectrum.export_caches()
    spectrum._default_spectrum(kvp=71.0, filtration_mm_al=2.5, anode_angle_deg=12.0)
    assert (71.0, 2.5, 12.0) not in spec_cache
