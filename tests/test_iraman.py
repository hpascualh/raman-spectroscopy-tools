from argparse import Namespace
from pathlib import Path
import numpy as np
import pytest

from plot_iraman_from_folder import (
    align_spectra,
    comparison_conditions,
    load_spectrum,
    preprocess,
    reject_z_outliers,
    shared_sample_y_limits,
)

ROOT = Path(__file__).resolve().parents[1]
SAMPLE = ROOT / 'examples' / 'iraman' / 'synthetic_day' / 'SYN_A_1EOD'


def args(**overrides):
    settings = dict(baseline=False, lam=1e5, p=1e-3, smooth=False,
                    norm='si', x_min=400, x_max=1800, no_std=False)
    settings.update(overrides)
    return Namespace(**settings)


def test_portable_file_parses_semicolon_shift_and_intensity():
    x, y = load_spectrum(SAMPLE / 'SP_1.txt')
    assert x.shape == y.shape
    assert x.size > 400
    assert np.all(np.diff(x) > 0)
    assert np.isfinite(y).all()


def test_duplicate_raman_shift_rejected(tmp_path):
    p = tmp_path / 'SP_1.txt'
    p.write_text('x;y\n510;1\n510;2\n520;3\n')
    with pytest.raises(ValueError, match='must be unique'):
        load_spectrum(p)


def test_alignment_interpolates_files_with_different_grids():
    paths = [SAMPLE / f'SP_{i}.txt' for i in range(1, 4)]
    x, y = align_spectra(paths)
    assert x.ndim == 1 and y.shape == (3, len(x))
    assert 200 < x[0] < x[-1] < 2000


def test_i_raman_preprocessing_normalizes_to_silicon_peak():
    x, values = align_spectra([SAMPLE / f'SP_{i}.txt' for i in range(1,4)])
    result = preprocess(x, values, args())
    np.testing.assert_allclose(result[:, (x >= 510) & (x <= 530)].max(axis=1), 1, atol=1e-12)


def test_outlier_filter_off_and_on():
    data = np.zeros((100, 3), dtype=float)
    data[-1, 1] = 42
    unchanged, removed = reject_z_outliers(data, None)
    assert unchanged.shape == (100, 3) and removed == 0
    cleaned, removed = reject_z_outliers(data, 8.0)
    assert cleaned.shape == (99, 3) and removed == 1


def test_common_y_limits_exclude_hidden_spectral_peaks():
    x = np.array([200, 400, 800, 1200, 1800, 2000], dtype=float)
    sample1 = ('a', x, np.array([10000, 0, 1, 2, 1, 10000]), np.full(6, .2))
    sample2 = ('b', x, np.array([10000, .5, 1, 3, 2, 10000]), np.full(6, .2))
    lo, hi = shared_sample_y_limits([sample1,sample2], args())
    assert lo <= -0.2 and 3.2 < hi < 10


def test_shared_condition_labels():
    common, labels = comparison_conditions(['SYN_A_1EOD_sensor','SYN_B_1EOD_sensor'])
    assert common == 'SYN 1EOD sensor'
    assert labels == ['A','B']
