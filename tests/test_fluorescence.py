"""蛍光X線（K殻のみ）のデータ層・サンプリング層・輸送組み込みのテスト。

docs/plan_fluorescence.md 参照。K殻蛍光のみ実装（L殻・カスケードは非目標）。
"""
from __future__ import annotations

import numpy as np

from chatcarlo.geometry import Geometry
from chatcarlo.materials import fluorescence_k_data, photo_element_weights
from chatcarlo.physics import _FLUOR_CUTOFF_KEV, sample_fluorescence
from chatcarlo.transport import transport_photons


def test_lead_k_edge_and_kalpha1():
    edge_keV, omega_k, line_energies, line_probs = fluorescence_k_data(82)
    assert abs(edge_keV - 88.00) < 0.1
    assert abs(omega_k - 0.96) < 0.02
    # KL3 = Kalpha1
    kalpha1 = line_energies[np.argmax(line_probs)]
    assert abs(kalpha1 - 74.97) < 0.1


def test_tungsten_k_edge_and_kalpha1():
    edge_keV, omega_k, line_energies, line_probs = fluorescence_k_data(74)
    assert abs(edge_keV - 69.53) < 0.1
    kalpha1 = line_energies[np.argmax(line_probs)]
    assert abs(kalpha1 - 59.32) < 0.1


def test_line_probs_sum_to_one():
    for z in (20, 26, 29, 74, 82):
        _, omega_k, line_energies, line_probs = fluorescence_k_data(z)
        if line_energies.size > 0:
            assert np.isclose(line_probs.sum(), 1.0)


def test_photo_element_weights_water_dominated_by_oxygen():
    zs, w = photo_element_weights("water", np.array([50.0]))
    assert np.allclose(w.sum(axis=0), 1.0)
    oxygen_idx = list(zs).index(8)
    assert w[oxygen_idx, 0] > 0.9


def test_no_fluorescence_below_k_edge():
    """K端未満のエネルギーでは蛍光は発生しない。"""
    n = 5000
    rng = np.random.default_rng(0)
    materials = np.full(n, "lead", dtype=object)
    e = np.full(n, 87.0)  # Pb K端 88.0 keV未満
    emit, e_line = sample_fluorescence(materials, e, rng)
    assert not np.any(emit)


def test_fluorescence_occurs_above_k_edge():
    n = 20_000
    rng = np.random.default_rng(1)
    materials = np.full(n, "lead", dtype=object)
    e = np.full(n, 100.0)  # Pb K端(88.0 keV)超
    emit, e_line = sample_fluorescence(materials, e, rng)
    assert np.mean(emit) > 0.1  # K殻分率×omega_kでそれなりの頻度のはず
    assert np.all(e_line[emit] >= _FLUOR_CUTOFF_KEV)
    assert np.all(e_line[emit] < 100.0)


def _thin_lead_slab_scene_arrays(n, energy_keV, seed):
    """薄い鉛スラブに垂直入射する単色光子束をtransport_photons直叩き用に組み立てる。"""
    thickness_cm = 0.5
    geometry = Geometry([{
        "name": "slab", "shape": "box", "material": "lead",
        "center": [0.0, 0.0, 0.0],
        "size_cm": [thickness_cm, 50.0, 50.0],
    }])
    rng = np.random.default_rng(seed)
    pos = np.tile(np.array([-thickness_cm / 2 - 5.0, 0.0, 0.0]), (n, 1))
    dirv = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    energy = np.full(n, float(energy_keV))
    return pos, dirv, energy, geometry, rng


def test_energy_conservation_with_fluorescence():
    n = 5000
    pos, dirv, energy, geometry, rng = _thin_lead_slab_scene_arrays(n, 100.0, seed=42)
    e_in_total = float(np.sum(energy))
    result = transport_photons(pos, dirv, energy, geometry, rng, fluorescence_enabled=True)
    e_deposited = sum(result.energy_deposited.values())
    e_escaped = float(np.sum(result.final_energy[result.escaped]))
    assert abs((e_deposited + e_escaped) - e_in_total) / e_in_total < 1e-9


def test_energy_conservation_without_fluorescence():
    n = 5000
    pos, dirv, energy, geometry, rng = _thin_lead_slab_scene_arrays(n, 100.0, seed=42)
    e_in_total = float(np.sum(energy))
    result = transport_photons(pos, dirv, energy, geometry, rng, fluorescence_enabled=False)
    e_deposited = sum(result.energy_deposited.values())
    e_escaped = float(np.sum(result.final_energy[result.escaped]))
    assert abs((e_deposited + e_escaped) - e_in_total) / e_in_total < 1e-9


def test_fluorescence_events_recorded_above_k_edge():
    n = 20_000
    pos, dirv, energy, geometry, rng = _thin_lead_slab_scene_arrays(n, 100.0, seed=7)
    result = transport_photons(pos, dirv, energy, geometry, rng, fluorescence_enabled=True)
    assert result.n_fluorescence > 0


def test_no_fluorescence_events_below_k_edge():
    n = 5000
    pos, dirv, energy, geometry, rng = _thin_lead_slab_scene_arrays(n, 87.0, seed=8)
    result = transport_photons(pos, dirv, energy, geometry, rng, fluorescence_enabled=True)
    assert result.n_fluorescence == 0


def test_fluorescence_disabled_produces_no_events():
    n = 20_000
    pos, dirv, energy, geometry, rng = _thin_lead_slab_scene_arrays(n, 100.0, seed=9)
    result = transport_photons(pos, dirv, energy, geometry, rng, fluorescence_enabled=False)
    assert result.n_fluorescence == 0
