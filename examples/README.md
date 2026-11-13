# Synthetic demonstration data

These files were randomly generated using fixed parameters and a deterministic seed. They are **not** measurements, do not represent bacterial signatures, and carry no experimental conclusions. They exist only to verify input parsing, preprocessing, plotting, and CLI examples.

- `xplora/SYN-Q001_.../Spot1.txt`, etc.: whitespace-delimited XploRA-style maps. First row: 2 coordinate placeholders plus 631 Raman-shift values. Remaining rows: 2 coordinate values + 631 spectral intensities.
- `iraman/synthetic_day/SYN_A_1EOD/SP_1.txt`, etc.: semicolon-separated i-Raman style TXT files with one header row, slight input-axis variations, and one spectrum per file.

Generated plots and summary CSVs appear in `results/` directories, which are ignored by Git.
