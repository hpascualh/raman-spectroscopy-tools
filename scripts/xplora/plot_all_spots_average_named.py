#!/usr/bin/env python3
"""Combine the spectra from every Spot map in one sample folder and plot mean +/- SD.

Keep this script in the same directory as plot_map_average_v2.py.
Example:
    python scripts/xplora/plot_all_spots_average_named.py "path/to/sample_folder"
"""

import argparse
import re
from pathlib import Path

import numpy as np
from plot_map_average_v2 import (load_xplora_map2array, remove_outlier_spectra,
                                 make_plot_title, plot_map_average)


def main():
    parser = argparse.ArgumentParser(description="Plot one average over all Spot maps in a sample folder")
    parser.add_argument("folder", type=Path, nargs="?", default=Path("."), help="Sample folder (default: current folder)")
    parser.add_argument("--z-thresh", type=float, default=8.0, help="Cosmic-ray outlier threshold (default: 8)")
    parser.add_argument("--no-outlier-filter", action="store_true")
    parser.add_argument("--smooth", action="store_true")
    parser.add_argument("--no-baseline", action="store_true")
    parser.add_argument("--no-norm-si", action="store_true")
    parser.add_argument("--no-std", action="store_true")
    parser.add_argument("--lam", type=float, default=10)
    parser.add_argument("--p", type=float, default=1e-4)
    args = parser.parse_args()

    folder = args.folder.resolve()
    if not folder.is_dir():
        parser.error(f"Folder not found: {folder}")
    maps = sorted((p for p in folder.glob("*.txt") if re.search(r"Spot\d+", p.stem, re.I)),
                  key=lambda p: p.name.lower())
    if not maps:
        parser.error(f"No Spot*.txt map files found in {folder}")

    all_spectra, common_shift = [], None
    total_removed = 0
    for path in maps:
        try:
            shift, spectra = load_xplora_map2array(
                path, smooth=args.smooth, baseline=not args.no_baseline,
                norm_si=not args.no_norm_si, LAM=args.lam, P=args.p)
            if common_shift is not None and (shift.shape != common_shift.shape or
                                             not np.allclose(shift, common_shift, rtol=0, atol=0.05)):
                raise ValueError(f"Raman-shift axis differs in {path.name}; cannot combine maps")
            if not args.no_outlier_filter:
                cleaned, _ = remove_outlier_spectra(spectra, z_thresh=args.z_thresh)
            else:
                cleaned = spectra
            removed = len(spectra) - len(cleaned)
            total_removed += removed
            print(f"{path.name}: kept {len(cleaned)}/{len(spectra)} spectra, removed {removed}")
            all_spectra.append(cleaned)
            if common_shift is None:
                common_shift = shift
        except (OSError, ValueError) as exc:
            parser.error(f"Problem processing {path.name}: {exc}")

    pooled = np.vstack(all_spectra)
    mean, std = pooled.mean(axis=0), pooled.std(axis=0)
    if not args.no_norm_si:
        si_mask = (common_shift >= 510) & (common_shift <= 530)
        si_max = mean[si_mask].max()
        if si_max <= 0:
            parser.error("Pooled mean has nonpositive Si peak")
        mean, std = mean / si_max, std / si_max

    title = re.sub(r"\s*--\s*Spot\s*\d+\s*$", "", make_plot_title(maps[0])) + " -- All spots"
    output_dir = Path.cwd() / "results" / "all_spots_average"
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_id = folder.parent.name + "_" + folder.name if re.search(r"[A-Za-z]{2,4}[-_]?Q\d{3,}", folder.parent.name, re.I) else folder.name
    output = output_dir / f"{sample_id}_all_spots_average.png"
    plot_map_average(common_shift, mean, std, title, output, show_std=not args.no_std)
    csv_output = output_dir / f"{sample_id}_all_spots_average.csv"
    np.savetxt(csv_output, np.column_stack((common_shift, mean, std)), delimiter=",",
               header="raman_shift,mean,std", comments="")
    print(f"Combined {len(maps)} maps, {len(pooled)} spectra retained, {total_removed} excluded")
    print(f"Saved: {output}")
    print(f"Saved: {csv_output}")


if __name__ == "__main__":
    main()
