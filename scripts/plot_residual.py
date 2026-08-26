#!/usr/bin/env python3
"""Plot residual convergence and force coefficients from a solver residue file.

Reads the Tecplot ASCII residue file written by ns_solver and produces:
  Left panel:  L2 residual norm vs iteration (semilogy)
  Right panel: Cl and Cd vs iteration (linear, twin axes)

Multiple files can be overlaid for comparing runs.

Usage:
    python3 scripts/plot_residual.py <residue.dat> [second.dat ...]
                                     [--labels "run1,run2"] [-o out.png]
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt


def read_residue(path):
    """Parse a residue .dat file, skipping header lines."""
    with open(path) as f:
        lines = f.readlines()
    data_lines = [
        ln.strip() for ln in lines
        if ln.strip() and not ln.strip().startswith(("TITLE", "VARIABLES"))
    ]
    data = np.array(
        [[float(x) for x in ln.split()] for ln in data_lines],
        dtype=np.float64,
    )
    return data


def main():
    ap = argparse.ArgumentParser(description="Plot residual convergence.")
    ap.add_argument("files", nargs="+", help="Residue .dat file(s)")
    ap.add_argument("--labels", default=None,
                     help="Comma-separated labels for each file")
    ap.add_argument("--start", type=int, default=1,
                    help="Skip first N iterations (default 1)")
    ap.add_argument("-o", "--output", default=None,
                    help="Output PNG path")
    args = ap.parse_args()

    if len(args.files) == 0:
        sys.exit("No input files specified")

    labels = args.labels.split(",") if args.labels else \
        [os.path.splitext(os.path.basename(f))[0] for f in args.files]
    if len(labels) != len(args.files):
        labels = [f"run {i+1}" for i in range(len(args.files))]

    fig, (ax_res, ax_force) = plt.subplots(
        1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [3, 2]}
    )

    colors = plt.cm.tab10(np.linspace(0, 0.9, len(args.files)))

    for idx, (path, label, color) in enumerate(zip(args.files, labels, colors)):
        if not os.path.isfile(path):
            sys.exit(f"File not found: {path}")

        data = read_residue(path)
        iters = data[:, 0].astype(int)
        res = data[:, 3]
        cl = data[:, 6]
        cd = data[:, 7]

        # Apply start offset
        mask = iters >= args.start
        it = iters[mask]
        r = res[mask]
        c_l = cl[mask]
        c_d = cd[mask]

        ax_res.semilogy(it, r, color=color, linewidth=1.2, label=label)
        ax_force.plot(it, c_l, color=color, linewidth=1.0, linestyle="-",
                      label=f"Cl ({label})" if len(args.files) == 1 else None)
        ax_force.plot(it, c_d, color=color, linewidth=1.0, linestyle="--",
                      label=f"Cd ({label})" if len(args.files) == 1 else None)

    ax_res.set_xlabel("Iteration")
    ax_res.set_ylabel(r"$\|\Delta \rho\|_2$")
    ax_res.set_title("Residual Convergence")
    ax_res.grid(True, alpha=0.4)
    ax_res.legend(fontsize=8)
    ax_res.set_xlim(left=args.start)

    ax_force.set_xlabel("Iteration")
    ax_force.set_ylabel("Coefficient")
    ax_force.set_title("Force Coefficients")
    ax_force.grid(True, alpha=0.4)
    if len(args.files) == 1:
        ax_force.legend(fontsize=8)

    base = os.path.splitext(os.path.basename(args.files[0]))[0]
    fig.suptitle(base.replace("_IsoCont", "").rsplit("Iter", 1)[0], fontsize=12)
    fig.tight_layout()

    out = args.output or f"{base}_residual.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
