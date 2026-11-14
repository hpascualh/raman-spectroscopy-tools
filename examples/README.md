# Synthetic demonstration data

These files were generated with a fixed random seed (20261113), varying Raman peak amplitudes so the plots show visible ±1 SD bands. They are **not** measurements, do not represent bacterial signatures, and carry no experimental conclusions. They exist only to verify input parsing, preprocessing, plotting, and CLI examples.

- `xplora/SYN-Q001_.../Spot1.txt`, etc.: whitespace-delimited XploRA-style maps. First row: 2 coordinate placeholders plus 631 Raman-shift values. Remaining rows: 2 coordinate values + 631 spectral intensities.
- `iraman/synthetic_day/SYN_A_1EOD/SP_1.txt`, etc.: semicolon-separated i-Raman style TXT files with one header row, slight input-axis variations, and one spectrum per file (six files in each of two sample folders).

Generated plots and summary CSVs appear in `results/` directories, which are ignored by Git.
