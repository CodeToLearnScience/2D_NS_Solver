#!/usr/bin/env python3
"""Plot contour fields from a solver solution (.dat) file.

Reads the legacy-compatible Tecplot ASCII format written by ns_solver and
produces filled-contour plots for density, pressure, Mach number, and
velocity magnitude.

Usage:
    python3 scripts/plot_contours.py <solution.dat> [--var mach] [-o out.png]

If --var is omitted, a 2×2 panel (rho, p, Mach, |V|) is produced.
Output PNG is saved next to the input file by default.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt


def read_solution(path):
    """Parse a Tecplot ASCII solution file into coordinates + field dict."""
    with open(path) as f:
        lines = f.readlines()

    # Parse header
    title_line = lines[0].strip()
    var_line = lines[1].strip()
    zone_line = lines[2].strip()

    # Extract variable names from "VARIABLES = "xc","yc","rho",..."
    vars_str = var_line.split("=", 1)[1]
    var_names = [v.strip().strip('"') for v in vars_str.split(",")]

    # Extract I and J from Zone line: 'Zone T="...", I=241, J=81'
    ni = nj = None
    for token in zone_line.replace(",", " ").split():
        if token.upper().startswith("I="):
            ni = int(token[2:])
        elif token.upper().startswith("J="):
            nj = int(token[2:])
    if ni is None or nj is None:
        raise ValueError(f"Cannot parse I/J from zone line: {zone_line}")

    # Read data (skip first 3 header lines)
    data_lines = lines[3:]
    data = np.array(
        [[float(x) for x in ln.split()] for ln in data_lines if ln.strip()],
        dtype=np.float64,
    )

    expected = ni * nj
    if data.shape[0] != expected:
        raise ValueError(f"Expected {expected} data points, got {data.shape[0]}")

    # Build dict: each column → 2D array (nj, ni), j-outer / i-inner
    fields = {}
    for col, name in enumerate(var_names):
        fields[name.lower()] = data[:, col].reshape(nj, ni)

    return fields, ni, nj, title_line


def main():
    ap = argparse.ArgumentParser(description="Plot solver solution contours.")
    ap.add_argument("file", help="Solution .dat file from ns_solver")
    ap.add_argument("--var", default=None,
                     choices=["rho", "p", "mach", "u", "v", "t"],
                     help="Plot single variable fullsize instead of 2×2 grid")
    ap.add_argument("-o", "--output", default=None,
                    help="Output PNG path (default: alongside input)")
    args = ap.parse_args()

    if not os.path.isfile(args.file):
        sys.exit(f"File not found: {args.file}")

    fields, ni, nj, title = read_solution(args.file)

    # Compute derived quantities
    rho = fields.get("rho")
    u = fields.get("u", np.zeros_like(rho))
    v = fields.get("v", np.zeros_like(rho))
    mach = fields.get("mach")
    if mach is None and rho is not None:
        gamma = 1.4
        p = fields.get("p", np.ones_like(rho))
        a = np.sqrt(gamma * p / rho)
        mach = np.sqrt(u * u + v * v) / a
    vel_mag = np.sqrt(u * u + v * v)

    x = fields.get("xc")
    y = fields.get("yc")

    base = os.path.splitext(os.path.basename(args.file))[0]

    if args.var:
        # Single-variable fullsize plot
        key = args.var.lower()
        data_map = {
            "rho": ("Density", rho),
            "p": ("Pressure", fields.get("p")),
            "mach": ("Mach Number", mach),
            "u": ("U-Velocity", u),
            "v": ("V-Velocity", v),
            "t": ("Temperature", fields.get("t")),
        }
        label, dat = data_map[key]
        fig, ax = plt.subplots(figsize=(10, 6))
        cf = ax.contourf(x, y, dat, 40, cmap="jet")
        fig.colorbar(cf, ax=ax, label=label)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(f"{label} — {base}")
        ax.set_aspect("equal")
        out = args.output or f"{base}_{key}.png"
    else:
        # Default 2×2 panel
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(base, fontsize=14)
        panels = [
            ("Density", rho, "jet"),
            ("Pressure [Pa]", fields.get("p"), "jet"),
            ("Mach Number", mach, "jet"),
            ("|Velocity| [m/s]", vel_mag, "jet"),
        ]
        for ax, (label, dat, cmap) in zip(axes.flat, panels):
            cf = ax.contourf(x, y, dat, 40, cmap=cmap)
            fig.colorbar(cf, ax=ax, label=label)
            ax.set_title(label)
            ax.set_aspect("equal")
        axes[0, 0].set_ylabel("y")
        axes[1, 0].set_ylabel("y")
        axes[1, 0].set_xlabel("x")
        axes[1, 1].set_xlabel("x")
        fig.tight_layout()
        out = args.output or f"{base}_contours.png"

    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
