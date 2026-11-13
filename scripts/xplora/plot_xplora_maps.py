#!/usr/bin/env python3
"""Plot the mean Raman spectrum (+/- 1 SD) from XploRA .txt map files."""

import argparse
import glob
import re
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from scipy.signal import savgol_filter
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


def correct_baseline(data, lam=10, p=1e-4, niter=10):
    """Asymmetric least-squares baseline, applied separately to each spectrum."""
    n = data.shape[1]
    D = diags([1., -2., 1.], [0, 1, 2], shape=(n-2, n), format="csc")
    penalty = lam * (D.T @ D)
    corrected = np.empty_like(data)
    for i, y in enumerate(data):
        w = np.ones(n)
        for _ in range(niter):
            z = spsolve(diags(w, format="csc") + penalty, w * y)
            w = np.where(y > z, p, 1-p)
        corrected[i] = y - z
    return corrected


def load_xplora_map2array(file_data, smooth=False, baseline=True, norm_si=True,
                          m=631, LAM=10, P=1e-4, window_length=11, polyorder=3,
                          si_range=(510, 530)):
    """Read XploRA data: first row = shifts; first two columns = map coordinates."""
    raman_shift = np.loadtxt(file_data, max_rows=1)
    data = np.loadtxt(file_data, skiprows=1, ndmin=2)
    values = data[:, 2:]

    if values.shape[1] != m:
        raise ValueError(f"Expected {m} spectral features; found {values.shape[1]}.")
    if raman_shift.size == data.shape[1]:
        raman_shift = raman_shift[2:]
    elif raman_shift.size != m:
        raise ValueError(f"Raman shifts ({raman_shift.size}) do not match intensities ({m}).")

    if baseline:
        values = correct_baseline(values, lam=LAM, p=P)
    if smooth:
        values = savgol_filter(values, window_length, polyorder, axis=1)
    if norm_si is not None and norm_si is not False:
        mask = (raman_shift >= si_range[0]) & (raman_shift <= si_range[1])
        if not mask.any():
            raise ValueError(f"No Raman-shift points within Si window {si_range} cm^-1.")
        si_peak = values[:, mask].max(axis=1)
        if np.any(si_peak <= 0):
            raise ValueError("At least one spectrum has a nonpositive Si peak.")
        values = values / si_peak[:, None]
    return raman_shift, values


def remove_outlier_spectra(spectra, z_thresh=8.0):
    """Remove spectra with unusually large pointwise z-scores (e.g., cosmic rays).

    For each wavenumber, compute z-scores across spectra. A spectrum is removed if
    any point exceeds |z| > z_thresh. The threshold is intentionally conservative.
    """
    if spectra.shape[0] < 3:
        keep_mask = np.ones(spectra.shape[0], dtype=bool)
        return spectra, keep_mask

    mean = spectra.mean(axis=0)
    std = spectra.std(axis=0, ddof=0)
    std_safe = np.where(std == 0, 1.0, std)
    z = (spectra - mean) / std_safe
    keep_mask = np.max(np.abs(z), axis=1) <= z_thresh

    # Avoid removing everything in pathological cases.
    if not keep_mask.any():
        keep_mask[:] = True
    return spectra[keep_mask], keep_mask


def make_plot_title(file_path):
    """Make a readable title using filename details, then parent-folder details."""
    path = Path(file_path)
    sources = (path.stem, path.parent.name)

    def find(pattern):
        for source in sources:
            match = re.search(pattern, source, flags=re.IGNORECASE)
            if match:
                return match
        return None

    sample = find(r"(?<![A-Za-z0-9])([A-Za-z]{2,4})[-_]?Q(\d{3,})(?![A-Za-z0-9])")
    eod = find(r"(?<![A-Za-z0-9])(\d+(?:\.\d+)?)\s*EOD(?![A-Za-z0-9])")
    conc = find(r"(?<![A-Za-z0-9])(\d+(?:\.\d+)?)\s*nM(?![A-Za-z0-9])")
    citrate = find(r"(?<![A-Za-z0-9])citrate(?![A-Za-z0-9])")
    sensor = find(r"(?<![A-Za-z0-9])LA[-_]?MUAM(?![A-Za-z0-9])")
    spot = find(r"(?<![A-Za-z0-9])Spot[-_ ]?(\d+)(?![A-Za-z0-9])")

    if not any((sample, eod, conc, citrate, sensor, spot)):
        return re.sub(r"Spot[-_ ]?(\d+)", r"Spot \1", path.stem.replace("_", " "), flags=re.IGNORECASE)

    title = f"{sample.group(1).upper()}-Q{sample.group(2)}" if sample else path.stem.replace("_", " ")
    if eod:
        title += f" at {eod.group(1)} EOD"
    if conc:
        title += f" with {conc.group(1)} nM"
        title += " Citrate" if citrate else " AuNPs"
    elif citrate:
        title += " with Citrate"
    if sensor:
        title += " on LA-MUAM sensor"
    if spot:
        title += f" -- Spot {spot.group(1)}"
    return title


def make_output_name(file_path):
    """Use parent folder name plus SpotX so output can be traced back to the experiment."""
    path = Path(file_path)
    spot = re.search(r"(?<![A-Za-z0-9])Spot[-_ ]?(\d+)(?![A-Za-z0-9])", path.stem, flags=re.IGNORECASE)
    if not spot:
        spot = re.search(r"(?<![A-Za-z0-9])Spot[-_ ]?(\d+)(?![A-Za-z0-9])", path.parent.name, flags=re.IGNORECASE)
    spot_part = f"_Spot{spot.group(1)}" if spot else ""
    return f"{path.parent.name}{spot_part}_average.png"


def plot_map_average(x, mean, std, title, output, show_std=True):
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(x, mean, linewidth=1.5)
    if show_std:
        ax.fill_between(x, mean-std, mean+std, color="#808080", alpha=0.35, edgecolor="none")
    ax.set(xlabel=r"Raman Shift (cm$^{-1}$)", ylabel="Intensity (A.U.)")
    ax.set_title(textwrap.fill(title, width=43, break_long_words=False), fontsize=11)

    ax.set_xlim(380, 2020)
    ax.set_xticks(np.arange(400, 2001, 200))
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.tick_params(axis="x", which="both", direction="in", top=False, labeltop=False,
                   bottom=True, labelbottom=True, labelsize=12)
    ax.tick_params(axis="y", which="both", direction="in", left=True, labelleft=True,
                   right=False, labelright=False, labelsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(output, dpi=250, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot mean +/- SD of XploRA Raman maps")
    parser.add_argument("files", nargs="*", help="Map .txt files (default: all .txt files in the current folder)")
    parser.add_argument("--title", help="Override the plot title (single file only)")
    parser.add_argument("--smooth", action="store_true", help="Apply Savitzky-Golay smoothing")
    parser.add_argument("--no-baseline", action="store_true", help="Skip baseline correction")
    parser.add_argument("--no-norm-si", action="store_true", help="Skip silicon normalization")
    parser.add_argument("--si-min", type=float, default=510, help="Lower bound of Si peak window in cm^-1 (default: 510)")
    parser.add_argument("--si-max", type=float, default=530, help="Upper bound of Si peak window in cm^-1 (default: 530)")
    parser.add_argument("--no-std", action="store_true", help="Hide standard-deviation shading")
    parser.add_argument("--lam", type=float, default=10, help="ALS baseline lambda (default: 10)")
    parser.add_argument("--p", type=float, default=1e-4, help="ALS baseline asymmetry (default: 1e-4)")
    parser.add_argument("--z-thresh", type=float, default=8.0, help="Outlier rejection threshold based on max |z| across a spectrum (default: 8.0)")
    parser.add_argument("--no-outlier-filter", action="store_true", help="Skip z-score outlier removal")
    args = parser.parse_args()

    paths = [Path(f) for pattern in args.files for f in glob.glob(pattern)] if args.files else sorted(Path.cwd().glob("*.txt"))
    paths = [p for p in paths if p.is_file()]
    if not paths:
        parser.error("No .txt map files found.")
    if args.title and len(paths) != 1:
        parser.error("--title may only be used with one input file.")

    output_dir = Path.cwd() / "results"
    output_dir.mkdir(exist_ok=True)

    for file in paths:
        try:
            x, spectra = load_xplora_map2array(
                file, smooth=args.smooth, baseline=not args.no_baseline,
                norm_si=not args.no_norm_si, LAM=args.lam, P=args.p,
                si_range=(args.si_min, args.si_max))

            n_before = spectra.shape[0]
            if not args.no_outlier_filter:
                spectra, keep_mask = remove_outlier_spectra(spectra, z_thresh=args.z_thresh)
                n_after = spectra.shape[0]
                print(f"{file.name}: kept {n_after}/{n_before} spectra; removed {n_before - n_after} outlier(s) with z-thresh={args.z_thresh:g}")
            else:
                print(f"{file.name}: kept {n_before}/{n_before} spectra; outlier filter OFF")

            mean, std = spectra.mean(axis=0), spectra.std(axis=0)
            if not args.no_norm_si:
                si_mask = (x >= args.si_min) & (x <= args.si_max)
                mean_si = mean[si_mask].max()
                if mean_si <= 0:
                    raise ValueError("Mean Si peak is not positive.")
                mean, std = mean / mean_si, std / mean_si
                print(f"Si window: {args.si_min:g}-{args.si_max:g} cm^-1, mean peak normalized to 1")

            output = output_dir / make_output_name(file)
            plot_map_average(x, mean, std, args.title or make_plot_title(file), output, show_std=not args.no_std)
            print(f"Saved: {output}")
        except (OSError, ValueError) as e:
            print(f"Skipped {file.name}: {e}")


if __name__ == "__main__":
    main()
