# Methods and scientific caveats

This document records **what the uploaded scripts do**. It is not a recommendation to harmonize different instrument settings or a validated Raman protocol.

## XploRA processing

1. `load_xplora_map2array` reads the first row as Raman shifts, each subsequent row as two map coordinates plus intensities, and requires 631 intensity columns by default.
2. Per-spectrum **ALS baseline subtraction** is applied by default with `lam=10`, `p=1e-4`, ten iterations. Smoothing is optional and uses Savitzky–Golay window length 11, polynomial order 3.
3. Each spectrum is normalized by the **maximum** corrected intensity in 510–530 cm⁻¹, after baseline and optional smoothing. Nonpositive silicon reference peaks fail validation.
4. `remove_outlier_spectra` computes pointwise z-scores across a map; any spectrum exceeding absolute z-score 8 at any spectral feature is dropped. If this would drop **all** spectra, it instead retains all; maps with fewer than three spectra are unchanged. The pooled-sample script filters **each map separately**, then concatenates retained spectra.
5. The per-map and pooled-sample scripts calculate mean and population-standard-deviation (`ddof=0`) across retained individual spectra and then **normalize the mean and SD again** by the mean's maximum in the silicon window. This changes scale but is part of current behavior.
6. The pooled-sample script requires matching Raman-shift grids between maps (absolute tolerance 0.05 cm⁻¹). It does not resample or interpolate map spectra.
7. The XploRA plots use a fixed 380–2020 cm⁻¹ range and shaded mean ± SD. Pooled CSV columns are `raman_shift,mean,std` and represent derived **sample averages**, not individual training spectra.
8. `compare_xplora_samples.py` requires directly matching saved shift grids. It offsets mean curves by `i * offset`, displays optional ±SD shading, and labels curves using the sample-name tokens that differ between files (or the full label when no tokens differ). Labels may be hidden with `--no-labels`.

## Portable i-Raman processing

1. Only files named `SP_<number>.txt` (case-insensitive) inside immediate subfolders of one experiment-day directory are imported. Files are expected to have a header followed by semicolon-delimited `Raman shift; intensity` records.
2. Input shift values are sorted, must be unique and finite, and every file must have at least three data rows.
3. The script finds the intersection of all supplied shift ranges **within each sample** and selects the densest input shift grid inside that intersection. Other spectra are linearly interpolated onto that grid.
4. Default ALS baseline parameters are `lam=1e5`, `p=1e-3`, ten iterations; optional Savitzky–Golay smoothing uses window 11 and polynomial order 3. By default the corrected spectrum is divided by its maximum in the 510–530 cm⁻¹ silicon window. Alternatives are max in displayed range or no normalization.
5. Outlier rejection is **off** unless `--z-thresh` is provided. Unlike the XploRA implementation, this routine raises an error if all spectra are rejected.
6. Mean and population SD are computed per sample, then optionally re-normalized by the sample mean's reference peak. Each CSV contains `raman_shift,mean,std`.
7. Individual plots use the **same y-axis limits**, calculated from visible 400–1800 cm⁻¹ mean ±SD values across samples; this avoids autoscale influence from hidden spectral regions. The day-level offset plot separately interpolates the sample means and SDs onto a common grid and labels each sample with its differing folder-name tokens.
8. The day-folder `results` location is determined from the **input path**, unlike the XploRA scripts, which write under the **current working directory**.

## Known limitations and interpretation

- The two scripts' ALS routines are mathematically similar but use **different defaults** and have not been refactored into a single implementation.
- Baseline correction and normalization parameters may need adjustment for other sample matrices, substrates, or Raman instruments. These defaults are carried forward for compatibility, not asserted as universally optimal.
- Standard deviation around pooled individual spectra should not be interpreted as uncertainty across independently prepared biological replicates.
- Interpolation produces comparable plotting grids but does **not** guarantee equal instrument resolution, sampling, response functions, or quantitative intensity comparability.
- Unsupervised analysis and supervised ML experiments are intentionally **outside** this repository; their validation and derived claims require separate treatment.
