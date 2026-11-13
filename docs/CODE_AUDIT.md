# Initial audit and migration decisions — 2026-11-13

## Source inventory

The four user-provided files were copied as the implementation authority; no numerical-processing or plotting logic was changed.

| Uploaded source | Repository path | Internal imports |
|---|---|---|
| `plot_map_average_v2(4).py` | `scripts/xplora/plot_map_average_v2.py` | Standard library, NumPy, Matplotlib, SciPy |
| `plot_all_spots_average_named(3).py` | `scripts/xplora/plot_all_spots_average_named.py` | Imports `load_xplora_map2array`, `remove_outlier_spectra`, `make_plot_title`, and `plot_map_average` from `plot_map_average_v2.py`; NumPy |
| `compare_sample_averages_from_csv(1).py` | `scripts/xplora/compare_sample_averages_from_csv.py` | Standard library, NumPy, Matplotlib |
| `plot_iraman_from_folder.py` | `scripts/iraman/plot_iraman_from_folder.py` | Standard library, NumPy, Matplotlib, SciPy |

The i-Raman input loader differs from the XploRA map reader. There is no need to create a common reader or Python package at this stage.

## Changes in this starter version

1. Removed download-version suffixes `(1)`, `(3)`, `(4)` when placing files in the repository, allowing the existing XploRA import to work directly.
2. Updated a stale filename and path in the pooled-XploRA script's top-level usage docstring.
3. Updated the i-Raman script's top-level example commands to reference its current filename and repository path.
4. Corrected the i-Raman `--offset` help string to reflect the existing numeric default (`2.0`).
5. Added installation documentation, synthetic examples, .gitignore, tests, and a methods reference.

**No default parameter value, function body, argument name, data format, plot style, or scientific computation was intentionally changed.**

## Deferred, targeted improvements

- Optional `--output-dir` flags for better location control. Existing scripts use different output conventions.
- The XploRA per-map script defaults to every `.txt` file in the working directory; use explicit filenames to avoid accidentally loading unrelated TXT files.
- The XploRA loader hardcodes `m=631` by default and only supports whitespace-delimited numeric export.
- The pooled-XploRA script matches `Spot<number>` in its `.txt` filenames, while the per-map script does not apply this filter to explicit filenames.
- The XploRA comparison script's default title contains the original experiment's conditions; the legend is disabled, reducing self-contained interpretability.
- File parsing and smoothing should be tested against approved, representative exports from both instruments.
- Add end-to-end testing on a Windows workstation before describing the repository as fully Windows validated.

## Boundary with ML research

The bacterial-ID ML scripts (unsupervised analysis, Round 1 classification, Round 2A/2B cross-sensor evaluations) will belong to a **separate, research-portfolio-oriented repository**. Those scripts are intentionally not needed to install or use these Raman processing tools. No ML source files have been included or modified here.

## Publication gate

No GitHub remote has been created; no sharing license has been chosen. Confirm Ragan Lab/UCI code ownership, collaborators' contributions, and data-release permissions before publishing publicly.

## Regression check completed on synthetic examples

The original uploaded scripts and packaged copies were executed with identical default arguments on the same synthetic inputs, using the tested local Python environment. The following outputs matched **byte-for-byte** (SHA-256), including PNGs: two pooled-XploRA CSVs and PNGs, the XploRA comparison PNG, two i-Raman per-sample CSVs and PNGs, and the i-Raman comparison PNG. The `plot_map_average_v2.py` scientific and plotting logic is unchanged; its CLI also passes an end-to-end synthetic smoke test.

These checks strongly support **behavior preservation for the tested examples and environment**, but do not establish equivalence for all experimental files or alternate dependency versions.
