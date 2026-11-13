# Raman spectroscopy processing and plotting tools

Small, script-based tools for processing and visualizing XploRA Raman maps and portable i-Raman single-spectrum exports. Intended as a **lab-oriented starter repository**, with minimal modifications to the original working scripts. This repository does **not** contain bacterial-identification ML code or experimental datasets.

> **Status:** Working starter version validated with **synthetic** sample files on Python 3.13 (NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8). **Not yet validated on real instrument exports or Windows.** Confirm instrument format, preprocessing settings, and scientific suitability for each dataset.

## 1. Install

Requirements: Python **3.11 or newer** (3.13 smoke-tested) and `pip`.

From the repository root, on **Windows PowerShell**:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

On **macOS/Linux**:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

If PowerShell blocks the activation script, you can run `.\.venv\Scripts\python.exe` in place of `python` (including for the commands below) without changing the system execution policy.

**Run commands from the repository root**, unless otherwise indicated. Input paths may be absolute or relative. The script files are deliberately kept together within each instrument folder so their imports work without packaging them into an installable Python library.

## 2. What's included

| Script | Input | Output |
|---|---|---|
| `scripts/xplora/plot_map_average_v2.py` | One or more raw XploRA map `.txt` files | Mean ± SD PNG **per map** in `results/` |
| `scripts/xplora/plot_all_spots_average_named.py` | One XploRA sample folder with `Spot<number>.txt` files | **Pooled-spectrum** mean ± SD PNG and `raman_shift,mean,std` CSV in `results/all_spots_average/` |
| `scripts/xplora/compare_sample_averages_from_csv.py` | Saved sample-average CSV files | Vertically offset comparison PNG in `results/all_spots_average/` (unless `--output` specified) |
| `scripts/iraman/plot_iraman_from_folder.py` | One experiment-day directory with sample subfolders containing `SP_<number>.txt` | Per-sample PNG and CSV plus a multi-sample offset PNG inside the **input day folder's** `results/` |

All scripts require only **NumPy, SciPy, and Matplotlib** (plus the Python standard library). There are no ML dependencies.

## 3. Try the synthetic examples

The files in `examples/` are **synthetically generated**, not lab measurements and not suitable for scientific interpretation.

### XploRA: per-map plot

```bash
python scripts/xplora/plot_map_average_v2.py "examples/xplora/SYN-Q001_1EOD_Citrate_LA_MUAM/Spot1.txt"
```

Produces `results/SYN-Q001_1EOD_Citrate_LA_MUAM_Spot1_average.png`.

### XploRA: pool every map in a sample

```bash
python scripts/xplora/plot_all_spots_average_named.py "examples/xplora/SYN-Q001_1EOD_Citrate_LA_MUAM"
python scripts/xplora/plot_all_spots_average_named.py "examples/xplora/SYN-Q002_1EOD_Citrate_LA_MUAM"
```

The pooled CSVs will be under `results/all_spots_average/`. The script combines **retained individual spectra** from all maps; it does **not** average map averages with equal weight.

### XploRA: compare saved sample-average CSVs

```bash
python scripts/xplora/compare_sample_averages_from_csv.py "results/all_spots_average/*_all_spots_average.csv" --title "Synthetic sample comparison" --offset 1.0
```

**Important:** the `--title` default still contains experiment-specific wording from the original script. Pass `--title` for other datasets. The comparison script currently draws each curve but keeps the Matplotlib legend disabled, so species/sample names are **not labeled on the figure** by default; CSV order determines vertical order.

### i-Raman: process an entire experiment day

```bash
python scripts/iraman/plot_iraman_from_folder.py "examples/iraman/synthetic_day"
```

Produces per-sample average CSVs and PNGs in `examples/iraman/synthetic_day/results/sample_averages/`, plus the day-level offset PNG in `examples/iraman/synthetic_day/results/`. (These folders are ignored by Git.)

The i-Raman script accepts `--no-baseline`, `--norm none|si|max`, `--smooth`, `--z-thresh`, `--x-min`, `--x-max`, `--offset`, and other options. To see all options for any script:

```bash
python scripts/iraman/plot_iraman_from_folder.py --help
```

### PowerShell wildcard note

Pass wildcard patterns **in quotes** to `compare_sample_averages_from_csv.py` and to the XploRA per-map script. The scripts expand the patterns themselves, making their behavior more consistent across Windows and macOS/Linux.

## 4. Expected raw input layouts

```text
YOUR_XPLORA_DATA/
  Sample_A/
    Spot1.txt
    Spot2.txt
  Sample_B/
    Spot1.txt

YOUR_IRAMAN_DAY/
  Sample_A/
    SP_1.txt
    SP_2.txt
  Sample_B/
    SP_1.txt
    SP_2.txt
```

**XploRA maps:** numeric whitespace-delimited text. First row holds Raman-shift values (either alone or preceded by two coordinate placeholders); remaining rows start with two map coordinates and then **631 spectral intensities**. This implementation requires 631 features per spectrum. For single-map plotting, pass explicit map filenames rather than relying on the default, which reads **every `.txt` file in the working directory**. For pooling, the loader searches for `.txt` filenames containing `Spot<number>` in the immediate folder (not recursively).

**Portable i-Raman:** `SP_1.txt`, `SP_2.txt`, etc. within each immediate sample subfolder; semicolon-delimited numeric shift and intensity, **one header row**, and optionally a UTF-8 BOM. The script interpolates individual files onto a common grid **within each sample**, then interpolates sample summaries onto a shared grid for the combined figure.

**Do not substitute CSV averages for raw individual spectra in ML training.** The CSV summaries are intended for plotting and examination of sample variability.

## 5. Preprocessing defaults (do not assume they are interchangeable)

| Parameter | XploRA scripts | i-Raman script |
|---|---|---|
| ALS baseline | On; `lam=10`, `p=1e-4` | On; `lam=1e5`, `p=1e-3` |
| Savitzky–Golay | Off by default; window 11, polynomial order 3 | Off by default; window 11, polynomial order 3 |
| Per-spectrum normalization | Silicon peak, 510–530 cm⁻¹ | Silicon peak by default; `--norm max` or `none` optional |
| Spectrum outlier removal | On; pointwise max absolute z-score threshold 8 | Off; optionally enable `--z-thresh` |
| Additional mean normalization | Yes, to the mean's Si peak after aggregation | Yes, to the mean's selected reference peak after aggregation |
| Presentation range | 380–2020 cm⁻¹ | 400–1800 cm⁻¹ |
| Pooled mean SD | Population SD over all retained spectra | Population SD over all retained spectra per sample |

See [`docs/METHODS.md`](docs/METHODS.md) for behavior, assumptions, and caveats. **Default values come from the provided scripts, not from a cross-instrument validated preprocessing protocol.**

## 6. Tests

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

Tests cover XploRA input parsing, normalization, outlier rejection, i-Raman parsing/alignment and plotting ranges, and all four command-line workflows on synthetic data. They do not establish equivalence on real measurements, scientific performance, or Windows portability.

## 7. Data management and publishing

- **Keep measured spectra outside this repository**, e.g. in a lab-managed data drive. Avoid committing raw experimental exports, identifiable information, cached arrays, or generated results.
- `.gitignore` excludes common local and derived-data paths; **review `git status` before committing**. `.gitignore` cannot prevent a previously tracked or arbitrarily named sensitive file from being committed.
- Example spectra are entirely synthetic and intentionally small.
- **No license is included yet.** Confirm ownership, collaborator permissions, institutional rules, and the preferred license before making the repository public or authorizing reuse. A public GitHub repository without an explicit license does not itself grant blanket reuse rights.

Before pushing to GitHub, follow the [publishing checklist](docs/PUBLISHING_CHECKLIST.md).

## 8. Planned small improvements (not implemented yet)

- CLI `--output-dir` options and a neutral default comparison title with clear sample labels.
- Targeted handling for maps that use other spectral-feature counts or slightly different raw exports.
- Real-file regression checks with approved, non-sensitive test material.
- Optional CI smoke tests and verified Windows instructions.

These should be introduced separately and checked against existing outputs, not folded into a wholesale rewrite.
