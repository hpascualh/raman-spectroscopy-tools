from pathlib import Path
import numpy as np
import pytest

from plot_map_average_v2 import (
    load_xplora_map2array,
    make_output_name,
    make_plot_title,
    remove_outlier_spectra,
)

ROOT = Path(__file__).resolve().parents[1]
SAMPLE = ROOT / 'examples' / 'xplora' / 'SYN-Q001_1EOD_Citrate_LA_MUAM' / 'Spot1.txt'


def test_map_loading_without_preprocessing():
    x, y = load_xplora_map2array(SAMPLE, baseline=False, norm_si=False)
    assert x.shape == (631,)
    assert y.shape == (12, 631)
    assert np.all(np.diff(x) > 0)
    assert np.isfinite(y).all()


def test_si_normalization_per_spectrum():
    x, y = load_xplora_map2array(SAMPLE, baseline=False, norm_si=True)
    mask = (x >= 510) & (x <= 530)
    np.testing.assert_allclose(np.max(y[:, mask], axis=1), 1.0, atol=1e-12)


def test_baseline_and_smoothing_defaults_remain_callable():
    x, y = load_xplora_map2array(SAMPLE, baseline=True, smooth=True)
    assert x.shape == (631,)
    assert y.shape == (12, 631)
    assert np.isfinite(y).all()


def test_map_rejects_wrong_number_of_features(tmp_path):
    f = tmp_path / 'Spot1.txt'
    # Deliberately only four wavenumbers, with two leading coordinate fields.
    np.savetxt(f, np.vstack(([0, 0, 510, 520, 530, 540], [1, 2, 1, 2, 3, 4])))
    with pytest.raises(ValueError, match='Expected 631'):
        load_xplora_map2array(f, baseline=False, norm_si=False)


def test_outlier_filter_retains_normal_spectra_and_rejects_spike():
    spectra = np.zeros((100, 631), dtype=float)
    spectra[-1, 200] = 50.0
    cleaned, kept = remove_outlier_spectra(spectra, z_thresh=8)
    assert cleaned.shape == (99, 631)
    assert kept.sum() == 99
    assert not kept[-1]


def test_map_title_and_output_name():
    assert 'SYN-Q001' in make_plot_title(SAMPLE)
    assert 'Spot 1' in make_plot_title(SAMPLE)
    assert make_output_name(SAMPLE).endswith('_Spot1_average.png')
