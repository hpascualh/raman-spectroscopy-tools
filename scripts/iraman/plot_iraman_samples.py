#!/usr/bin/env python3
"""Process single-spectrum SP_<number>.txt Raman files from one experiment day.

From the repository root, for example:
    python scripts/iraman/plot_iraman_samples.py "20261105_experiment"
    python scripts/iraman/plot_iraman_samples.py "20261105_experiment" --no-baseline --norm max
    python scripts/iraman/plot_iraman_samples.py "20261105_experiment" --offset 0.8 --no-std

Each immediate sample subfolder is processed only if it contains SP_<number>.txt files.
No Excel conversion or ordering/continuous spectrum indices are required.
"""
import argparse
from collections import Counter
import re
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, MultipleLocator
from scipy.signal import savgol_filter
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

SP_PATTERN = re.compile(r"SP_\d+\.txt", flags=re.IGNORECASE)


def load_spectrum(path):
    """Read two semicolon-separated columns: Raman shift; dark-subtracted intensity."""
    data = np.loadtxt(path, delimiter=";", skiprows=1, encoding="utf-8-sig", ndmin=2)
    if data.ndim != 2 or data.shape[1] != 2 or len(data) < 3 or not np.isfinite(data).all():
        raise ValueError(f"Unexpected data or non-finite values in {path}")
    x, y = data.T
    order = np.argsort(x)
    x, y = x[order], y[order]
    if np.any(np.diff(x) <= 0):
        raise ValueError(f"Raman-shift values must be unique in {path}")
    return x, y


def align_spectra(files):
    """Align spectra on the densest supplied grid within their shared shift range."""
    series = [load_spectrum(f) for f in files]
    lo, hi = max(x[0] for x, _ in series), min(x[-1] for x, _ in series)
    if lo >= hi:
        raise ValueError("These SP files have no overlapping Raman-shift range.")
    ref_x = max((x for x, _ in series), key=len)
    common_x = ref_x[(ref_x >= lo) & (ref_x <= hi)]
    if common_x.size < 3:
        raise ValueError("Fewer than three Raman-shift points are shared by these files.")
    aligned = np.vstack([np.interp(common_x, x, y) for x, y in series])
    return common_x, aligned


def baseline_als(values, lam=1e5, p=1e-3, iterations=10):
    """Asymmetric least-squares baseline subtraction per spectrum."""
    n = values.shape[1]
    d = diags([1., -2., 1.], [0, 1, 2], shape=(n-2, n), format="csc")
    penalty = lam * (d.T @ d)
    result = np.empty_like(values)
    for i, y in enumerate(values):
        weights = np.ones(n)
        for _ in range(iterations):
            z = spsolve(diags(weights, format="csc") + penalty, weights * y)
            weights = np.where(y > z, p, 1-p)
        result[i] = y - z
    return result


def preprocess(x, spectra, args):
    if args.baseline:
        spectra = baseline_als(spectra, lam=args.lam, p=args.p)
    if args.smooth:
        spectra = savgol_filter(spectra, 11, 3, axis=1)
    if args.norm == "si":
        mask = (x >= 510) & (x <= 530)
        if not mask.any():
            raise ValueError("Si normalization requested but no points from 510 to 530 cm^-1 found.")
        scale = spectra[:, mask].max(axis=1)
    elif args.norm == "max":
        window = (x >= args.x_min) & (x <= args.x_max)
        if not window.any():
            raise ValueError("No Raman-shift points in selected plot range for max normalization.")
        scale = spectra[:, window].max(axis=1)
    else:
        return spectra
    if np.any(scale <= 0):
        raise ValueError("Cannot normalize a spectrum with a nonpositive reference peak.")
    return spectra / scale[:, None]


def reject_z_outliers(spectra, threshold):
    """Conservative old-style pointwise z-score filter; optional, not enabled by default."""
    if threshold is None or len(spectra) < 3:
        return spectra, 0
    mean, std = spectra.mean(axis=0), spectra.std(axis=0)
    z = np.abs((spectra - mean) / np.where(std > 0, std, 1))
    keep = np.max(z, axis=1) <= threshold
    if not keep.any():
        raise ValueError("Outlier filtering rejected all spectra. Use a higher --z-thresh.")
    return spectra[keep], int((~keep).sum())


def plot_axes(ax, x, args, *, normalized=False):
    ax.set_xlabel(r"Raman Shift (cm$^{-1}$)", fontsize=11)
    ax.set_ylabel("Intensity (A.U.)", fontsize=11)
    ax.set_xlim(args.x_min, args.x_max)
    ticks = np.arange(200, args.x_max + 201, 200)
    ticks = ticks[(ticks >= args.x_min) & (ticks <= args.x_max)]
    if len(ticks):
        ax.set_xticks(ticks)
        ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    if normalized:
        ax.yaxis.set_major_locator(MultipleLocator(1))
    else:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
    ax.tick_params(axis="both", which="both", direction="in", labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def create_sample_plot(x, mean, std, sample_name, out, args, y_limits=None):
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(x, mean, linewidth=1.5)
    if not args.no_std:
        ax.fill_between(x, mean-std, mean+std, color="#808080", alpha=0.35, linewidth=0)
    plot_axes(ax, x, args, normalized=args.norm != "none")
    if y_limits is not None:
        ax.set_ylim(*y_limits)
    ax.set_title(textwrap.fill(sample_name.replace("_", " "), width=45, break_long_words=False), fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=250, bbox_inches="tight")
    plt.close(fig)


def shared_sample_y_limits(samples, args):
    """Common limits from visible mean +/- std across all individual samples.

    Matplotlib otherwise uses all x-values when autoscaling y, even if xlim
    shows only a smaller spectral range. This intentionally excludes data
    outside [x_min, x_max] and does not apply to the offset comparison plot.
    """
    lowest, highest = np.inf, -np.inf
    for _, x, mean, std in samples:
        visible = (x >= args.x_min) & (x <= args.x_max)
        if not visible.any():
            raise ValueError("A sample has no Raman-shift points in the displayed range.")
        shading = std[visible] if not args.no_std else np.zeros(visible.sum())
        lowest = min(lowest, np.min(mean[visible] - shading))
        highest = max(highest, np.max(mean[visible] + shading))

    span = max(highest - lowest, 0.1)
    ymin = min(0.0, lowest - 0.05*span)
    ymax = highest + 0.04*span # 0.08 og
    if args.norm != "none":
        ymin = np.floor(ymin*2)/2 if ymin < 0 else 0.0
        ymax = np.ceil(ymax)  # tidy integer upper bound for normalized spectra
    if ymax <= ymin:
        ymax = ymin + 1
    return float(ymin), float(ymax)


def comparison_conditions(names):
    """Separate tokens shared by every sample folder from sample-specific tokens.

    Underscores and spaces mark condition boundaries; internal hyphens and
    decimal points are preserved (e.g. EC-Q001, 0.325nM, 3-5p).
    """
    parts = [re.split(r"[_\s]+", name.strip("_ ")) for name in names]
    counts = [Counter(token.casefold() for token in part if token) for part in parts]
    common = counts[0].copy()
    for count in counts[1:]:
        common &= count

    title_parts = []
    remaining = common.copy()
    for token in parts[0]:
        key = token.casefold()
        if remaining[key] > 0:
            title_parts.append(token)
            remaining[key] -= 1

    unique_labels = []
    for original, tokens in zip(names, parts):
        remaining = common.copy()
        differences = []
        for token in tokens:
            key = token.casefold()
            if remaining[key] > 0:
                remaining[key] -= 1
            else:
                differences.append(token)
        unique_labels.append(" ".join(differences) or original.replace("_", " "))
    return " ".join(title_parts), unique_labels


def create_comparison(samples, out, args):
    lo = max(x.min() for _, x, _, _ in samples)
    hi = min(x.max() for _, x, _, _ in samples)
    if lo >= hi:
        raise ValueError("Samples do not share a Raman-shift range for comparison.")
    ref_x = max((x for _, x, _, _ in samples), key=len)
    x = ref_x[(ref_x >= lo) & (ref_x <= hi)]
    x = x[(x >= args.x_min) & (x <= args.x_max)]
    if x.size < 3:
        raise ValueError("Not enough common Raman-shift points in the requested plot range.")

    means = [np.interp(x, sx, mean) for _, sx, mean, _ in samples]
    stds = [np.interp(x, sx, std) for _, sx, _, std in samples]
    offset = args.offset
    if offset <= 0:
        raise ValueError("--offset must be positive.")

    names = [name for name, _, _, _ in samples]
    common_conditions, labels = comparison_conditions(names)
    fig, ax = plt.subplots(figsize=(7, len(samples)))
    # By default labels sit just beyond the right edge of the plotted spectra;
    # a fixed point-based offset prevents overlap with the lines themselves.
    label_x = args.label_x if args.label_x is not None else float(x[-1])
    if not (x[0] <= label_x <= x[-1]):
        raise ValueError(f"--label-x must be between {x[0]:g} and {x[-1]:g} cm^-1")
    y_values = []
    for i, (label, mean, std) in enumerate(zip(labels, means, stds)):
        y = mean + i*offset
        y_values.append(y)
        ax.plot(x, y, linewidth=1.4)
        if not args.no_std:
            ax.fill_between(x, y-std, y+std, alpha=0.17, color="#808080", linewidth=0)
        label_at_x = float(np.interp(label_x, x, y))
        ax.annotate(textwrap.fill(label, width=32, break_long_words=False),
                    xy=(label_x, label_at_x), xytext=(9, args.label_gap),
                    textcoords="offset points", ha="left", va="bottom",
                    fontsize=8, annotation_clip=False)

    plot_axes(ax, x, args, normalized=args.norm != "none")
    ax.set_ylabel("Intensity (A.U.)", fontsize=11)
    title = args.title or common_conditions or args.day_folder.name.replace("_", " ")
    ax.set_title(textwrap.fill(title, width=64, break_long_words=False), fontsize=11)
    # Extend y range so the labels of the uppermost spectrum remain visible.
    ymin = min(np.min(y - s) for y, s in zip(y_values, stds))
    ymax = max(np.max(y + s) for y, s in zip(y_values, stds))
    ax.set_ylim(ymin - 0.06*offset, ymax + 0.17*offset)
    fig.tight_layout()
    fig.savefig(out, dpi=250, bbox_inches="tight")
    plt.close(fig)
    return offset


def main():
    parser = argparse.ArgumentParser(description="Average SP_i.txt Raman spectra for each sample and compare samples")
    parser.add_argument("day_folder", type=Path, help="Experiment-day directory containing sample subfolders")
    parser.add_argument("--offset", type=float, default=2.0, help="Vertical offset (default: 2)")
    parser.add_argument("--x-min", type=float, default=400, help="Left plot limit (default: 400 cm^-1)")
    parser.add_argument("--x-max", type=float, default=1800, help="Right plot limit (default: 1800 cm^-1)")
    parser.add_argument("--label-x", type=float, default=None, help="Raman shift where labels appear (default: near right edge)")
    parser.add_argument("--label-gap", type=float, default=8, help="Vertical label clearance above line endpoints, in points (default: 8)")
    parser.add_argument("--title", default=None, help="Combined-figure title")
    parser.add_argument("--no-std", action="store_true", help="Hide standard-deviation shading")
    baseline_group = parser.add_mutually_exclusive_group()
    baseline_group.add_argument("--baseline", dest="baseline", action="store_true", help="Use ALS baseline correction (default)")
    baseline_group.add_argument("--no-baseline", dest="baseline", action="store_false", help="Skip ALS baseline correction")
    parser.set_defaults(baseline=True)
    parser.add_argument("--lam", type=float, default=1e5, help="ALS lambda (default 1e5)")
    parser.add_argument("--p", type=float, default=1e-3, help="ALS p (default 1e-3)")
    parser.add_argument("--smooth", action="store_true", help="Apply Savitzky-Golay smoothing (off by default)")
    parser.add_argument("--norm", choices=("none", "si", "max"), default="si",
                        help="Normalize each spectrum by Si peak (default), maximum in plotted range, or not at all")
    parser.add_argument("--z-thresh", type=float, default=None,
                        help="Optional per-wavenumber z-score spectrum rejection threshold (e.g., 8)")
    args = parser.parse_args()
    day = args.day_folder.resolve()
    if not day.is_dir():
        parser.error(f"Experiment-day folder not found: {day}")
    if args.z_thresh is not None and args.z_thresh <= 0:
        parser.error("--z-thresh must be positive")
    if args.x_min >= args.x_max:
        parser.error("--x-min must be less than --x-max")
    if args.offset <= 0:
        parser.error("--offset must be positive")
    if args.label_gap < 0:
        parser.error("--label-gap cannot be negative")

    results = day / "results"
    individual_dir = results / "sample_averages"
    individual_dir.mkdir(parents=True, exist_ok=True)
    samples = []
    for folder in sorted(day.iterdir(), key=lambda f: f.name.lower()):
        if not folder.is_dir() or folder == results:
            continue
        files = sorted((f for f in folder.iterdir() if f.is_file() and SP_PATTERN.fullmatch(f.name)),
                       key=lambda f: int(f.stem.split("_")[1]))
        if not files:
            continue
        try:
            x, spectra = align_spectra(files)
            spectra = preprocess(x, spectra, args)
            spectra, removed = reject_z_outliers(spectra, args.z_thresh)
            mean, std = spectra.mean(axis=0), spectra.std(axis=0)
            if args.norm != "none":
                mask = ((x >= 510) & (x <= 530)) if args.norm == "si" else ((x >= args.x_min) & (x <= args.x_max))
                scale = mean[mask].max()
                if scale <= 0:
                    raise ValueError("Mean normalization reference peak is nonpositive.")
                mean, std = mean/scale, std/scale

            prefix = folder.name
            png_out = individual_dir / f"{prefix}_average.png"
            csv_out = individual_dir / f"{prefix}_average.csv"
            np.savetxt(csv_out, np.column_stack([x, mean, std]), delimiter=",",
                       header="raman_shift,mean,std", comments="")
            samples.append((prefix, x, mean, std))
            print(f"{folder.name}: {len(files)} spectra found, {len(spectra)} retained, "
                  f"{removed} excluded -> {csv_out.name}")
        except (OSError, ValueError) as e:
            print(f"WARNING: Skipped {folder.name}: {e}")

    if not samples:
        parser.error(f"No usable sample folders containing SP_<number>.txt in {day}")

    # Generate every individual plot with the SAME y-axis limits based only
    # on visible data (including +/- std unless --no-std was requested).
    y_limits = shared_sample_y_limits(samples, args)
    print(f"Shared individual-plot y-axis: {y_limits[0]:g} to {y_limits[1]:g} "
          f"(calculated within {args.x_min:g}-{args.x_max:g} cm^-1)")
    for name, x, mean, std in samples:
        png_out = individual_dir / f"{name}_average.png"
        create_sample_plot(x, mean, std, name, png_out, args, y_limits=y_limits)
        # create_sample_plot(x, mean, std, name, png_out, args, y_limits=None)
        print(f"Saved: {png_out}")

    if len(samples) >= 2:
        combined_out = results / f"{day.name}_comparison_offset.png"
        try:
            offset = create_comparison(samples, combined_out, args)
            print(f"Combined {len(samples)} samples with offset {offset:g} -> {combined_out}")
        except ValueError as e:
            print(f"WARNING: Could not create combined plot: {e}")
    else:
        print("One sample found, so no multi-sample comparison was generated.")
    print(f"Results saved under: {results}")


if __name__ == "__main__":
    main()
