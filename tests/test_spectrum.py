"""線源スペクトル生成（SpekPy統合）のテスト。"""
from __future__ import annotations

import numpy as np

from vivemonte import spectrum


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
