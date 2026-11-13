#!/usr/bin/env python3
"""Offset plot comparing processed sample averages from saved CSV files.

Examples
--------
python compare_sample_averages_from_csv.py "results/all_spots_average/*_all_spots_average.csv" --offset 0.8
python compare_sample_averages_from_csv.py \
    "results/all_spots_average/sample1_all_spots_average.csv" \
    "results/all_spots_average/sample2_all_spots_average.csv" \
    --offset 1.0 --title "4 bacterial species"
"""

import argparse
import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


def pretty_label(csv_path: Path) -> str:
    name = csv_path.stem
    if name.endswith("_all_spots_average"):
        name = name[:-len("_all_spots_average")]
    return name


def main():
    parser = argparse.ArgumentParser(description="Offset plot of pooled sample means from saved CSV files")
    parser.add_argument("csv_files", nargs="+", help="CSV files or wildcard patterns, in bottom-to-top plot order")
    parser.add_argument("--offset", type=float, default=0.8, help="Vertical distance between spectra (default: 0.8)")
    parser.add_argument("--no-std", action="store_true", help="Hide standard-deviation shading")
    parser.add_argument("--title", default="4 species comparison - 1 EOD, 0.325 nM", help="Figure title")
    parser.add_argument("--output", default=None, help="Output PNG path (default: results/all_spots_average/sample_comparison_offset.png)")
    args = parser.parse_args()

    paths = []
    for pattern in args.csv_files:
        matches = [Path(p) for p in glob.glob(pattern)]
        if matches:
            paths.extend(matches)
        else:
            p = Path(pattern)
            if p.is_file():
                paths.append(p)

    # deduplicate while preserving order
    unique_paths = []
    seen = set()
    for p in paths:
        rp = str(p.resolve())
        if rp not in seen:
            unique_paths.append(p)
            seen.add(rp)
    paths = unique_paths

    if not paths:
        parser.error("No CSV files found.")

    samples = []
    for csv_path in paths:
        if not csv_path.is_file():
            parser.error(f"Missing CSV file: {csv_path}")
        array = np.loadtxt(csv_path, delimiter=",", skiprows=1, ndmin=2)
        if array.shape[1] != 3:
            parser.error(f"Expected 3 columns (raman_shift,mean,std) in {csv_path}")
        x, mean, std = array.T
        samples.append((pretty_label(csv_path), x, mean, std))

    common_x = samples[0][1]
    for name, x, _, _ in samples[1:]:
        if x.shape != common_x.shape or not np.allclose(x, common_x, rtol=0, atol=0.05):
            parser.error(f"Raman shifts in {name} do not match the first sample; cannot overlay directly")

    fig, ax = plt.subplots(figsize=(6, 5))
    for i, (name, x, mean, std) in enumerate(samples):
        y = mean + i * args.offset
        ax.plot(x, y, linewidth=1.5, label=name)
        if not args.no_std:
            ax.fill_between(x, y - std, y + std, color="#808080", alpha=0.2, linewidth=0)

    ax.set(xlabel=r"Raman Shift (cm$^{-1}$)", ylabel="Normalized intensity + offset (A.U.)")
    ax.set_title(args.title, fontsize=11)
    ax.set_xlim(380, 2020)
    ax.set_xticks(np.arange(400, 2001, 200))
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.tick_params(axis="x", which="both", direction="in", top=False, labeltop=False,
                   bottom=True, labelbottom=True, labelsize=10)
    ax.tick_params(axis="y", which="both", direction="in", left=True, labelleft=True,
                   right=False, labelright=False, labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # ax.legend()
    # ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=9)
    fig.tight_layout()

    default_out = Path("results") / "all_spots_average" / "sample_comparison_offset.png"
    out = Path(args.output) if args.output else default_out
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
