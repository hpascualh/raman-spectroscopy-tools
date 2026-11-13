"""Exercise four original command-line entry points on synthetic instrument files."""
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts'
EXAMPLES = ROOT / 'examples'


def run_cli(cwd, script, *arguments):
    completed = subprocess.run(
        [sys.executable, str(SCRIPTS / script), *map(str, arguments)],
        cwd=cwd, text=True, capture_output=True, check=False,
    )
    assert completed.returncode == 0, (completed.stdout, completed.stderr)
    return completed.stdout


def test_all_four_cli_workflows(tmp_path):
    xplora = EXAMPLES / 'xplora'
    sample1 = xplora / 'SYN-Q001_1EOD_Citrate_LA_MUAM'
    sample2 = xplora / 'SYN-Q002_1EOD_Citrate_LA_MUAM'

    run_cli(tmp_path, 'xplora/plot_map_average_v2.py', sample1 / 'Spot1.txt')
    assert (tmp_path / 'results' / f'{sample1.name}_Spot1_average.png').is_file()

    for sample in (sample1, sample2):
        run_cli(tmp_path, 'xplora/plot_all_spots_average_named.py', sample)
        csv = tmp_path / 'results' / 'all_spots_average' / f'{sample.name}_all_spots_average.csv'
        arr = np.loadtxt(csv, delimiter=',', skiprows=1)
        assert arr.shape == (631, 3)
        assert np.isfinite(arr).all()
        assert csv.with_suffix('.png').is_file()

    run_cli(tmp_path, 'xplora/compare_sample_averages_from_csv.py',
            str(tmp_path / 'results' / 'all_spots_average' / '*_all_spots_average.csv'),
            '--title', 'Synthetic comparison')
    assert (tmp_path / 'results' / 'all_spots_average' / 'sample_comparison_offset.png').is_file()

    # i-Raman script writes inside its input folder, so copy demo data to temp.
    day = tmp_path / 'synthetic_day'
    shutil.copytree(EXAMPLES / 'iraman' / 'synthetic_day', day,
                    ignore=shutil.ignore_patterns('results'))
    run_cli(tmp_path, 'iraman/plot_iraman_from_folder.py', day)
    assert (day / 'results' / 'synthetic_day_comparison_offset.png').is_file()
    for name in ('SYN_A_1EOD', 'SYN_B_1EOD'):
        csv = day / 'results' / 'sample_averages' / f'{name}_average.csv'
        arr = np.loadtxt(csv, delimiter=',', skiprows=1)
        assert arr.shape[1] == 3
        assert np.isfinite(arr).all()
        assert csv.with_suffix('.png').is_file()
