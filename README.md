# 2D_NS_Solver

2D compressible Euler / Navier–Stokes finite-volume solver on body-fitted
structured grids. This repository is mid-refactor: the original single-binary
C++98 code (`src/*.cpp` legacy tree, target `ns_legacy` → binary
`2D_NS_Solver`) is being rebuilt as a modern C++23 stack (targets `ns_core`,
`ns_solver`, `ns_config`) with TOML configuration, Google Test coverage, and
optional MPI parallelization. See `REFACTORING_PLAN.md` for the full roadmap,
audit findings, and per-phase verification gates.

## Status

Phases 0–7 complete: audit & bug fixes, CMake+GTest scaffolding, TOML config
layer, ghost-aware fields/mesh/IO infrastructure, physics + inviscid scheme
ports with bit-for-bit legacy parity goldens, single-rank end-to-end driver,
and j-slab MPI decomposition.

Quick summary of what runs today:

| Target | What it is |
|---|---|
| `ns_legacy` / `2D_NS_Solver` | Frozen pre-refactor solver (parity reference) |
| `ns_core` | Fields, mesh/metrics, numerics, physics, IO, MPI library |
| `ns_config` | Typed TOML config loader + legacy `.cfg` importer |
| `ns_solver` | New driver: `ns_solver <config.toml> <grid.dat> <out-dir> [iters] [--restart file]` |
| `cfg2toml` | Converts legacy `Input*.cfg` + `GridTop*.dat` pairs to canonical TOML |

## Building

```bash
cmake --preset release          # or: debug / asan
cmake --build --preset release -j
ctest --preset release          # unit + parity tests
```

MPI flavour (j-slab domain decomposition):

```bash
cmake -S . -B build-mpi -D CMAKE_BUILD_TYPE=Release -D NS_ENABLE_MPI=ON
cmake --build build-mpi -j
mpirun -np 4 build-mpi/ns_solver configs/RampEuler.toml \
    input/RampGrid24080.dat out-dir 200
```

`NS_ENABLE_NATIVE_ARCH=ON` adds `-march=native`; note this enables FMA
contraction and shifts last-bit results — parity gates must run without it.

## Configuration

Canonical schema lives in `configs/*.toml` (migrated from `input/Input*.cfg`
by `scripts/migrate_configs.sh`). Boundary segments keep the legacy node-index
convention; edges are named `imin/imax/jmin/jmax`. Restart:

```bash
ns_solver cfg.toml grid.dat out 100            # writes nothing special
ns_solver cfg.toml grid.dat out 10 --restart checkpoint.rst   # continue
```

Checkpoints are written by adding `solver.write_restart(path)` calls (see
`src/solver/euler_solver.hpp`).

## Testing & parity

- `tests/` — 66 GTest cases: math utilities, config load/validation, legacy
  import round-trips, fields/ghost semantics, metrics vs hand-computed
  geometry, restart bit-exactness, kernel-level **bitwise** golden parity with
  the legacy solver (physics, LLF/Roe/NKFDS-MOVERS chains), MPI halo exchange.
- `scripts/capture_regression.sh` / `compare_regression.sh` — full-file
  regression baselines against the frozen legacy binary.
- `scripts/check_e2e_parity.sh` — serial new-vs-legacy trajectory check.
- `scripts/check_mpi_parity.sh [iters]` — np=1/2/4 vs legacy.
- `scripts/bench.sh [iters]` — timings: legacy vs new serial vs MPI.

Known numerical caveat (documented in `REFACTORING_PLAN.md` §6): the legacy
NKFDS-MOVERS scheme quantizes wave speeds to 3 decimals; sub-roundoff codegen
differences therefore amplify to ~1e-4 relative field deviations over long
runs. Kernel-level parity gates remain bitwise.

## Layout

```
src/config/    TOML loading, validation, legacy importer
src/core/      error type
src/fields/    Field<T>/MultiField<T>: contiguous, ghost-aware, View2D
src/mesh/      StructuredMesh, legacy grid reader (incl. slab loader), metrics
src/numerics/  face flux data, LLF, NKFDS-MOVERS, Roe schemes
src/physics/   IdealGas EOS, Sutherland transport, dependent variables,
               Green-Gauss face gradients, viscous flux
src/io/        VTK writer, binary restart, (legacy-format writers in solver)
src/parallel/  SlabDecomposition, HaloExchange
src/solver/    EulerSolver driver (serial + MPI)
tools/cfg2toml config migration tool
legacy refs:   src/*.cpp (frozen), input/*.cfg|*.dat (cases)
```

## Scripts

| Script | Purpose |
|---|---|
| `scripts/migrate_configs.sh` | Convert all convertible legacy inputs to TOML |
| `scripts/capture_regression.sh` | Capture golden outputs from legacy binary |
| `scripts/check_e2e_parity.sh` | Serial new-vs-legacy comparison |
| `scripts/check_mpi_parity.sh` | np=1/2/4 vs legacy residuals |
| `scripts/bench.sh` | Timing comparison across implementations |

## Validation test suite

Nine published benchmark cases from the two papers behind the numerical
methods — Maruthi & Rao, *An Entropy Stable Central Solver for Euler
Equations* (arXiv:1504.06442) and Shrinath et al., *A Kinetic Flux Difference
Splitting method for compressible flows* (Computers & Fluids 250, 2023) —
are configured under `configs/`, one TOML per case. Run any of them with

```bash
./build/release/ns_solver configs/<Case>.toml input/<grid>.dat output/<dir> [max-iters]
python3 scripts/plot_contours.py  output/<dir>/*IsoCont*.dat   # contour plots
python3 scripts/plot_residual.py  output/<dir>/Residue*.dat    # convergence
```

Reference plots produced with these exact commands are stored in each case's
`output/<case>/` directory. All Euler cases are nondimensional with
u* = M∞ and p* = 1/γ; NS cases use the legacy NS scaling (T*=u*=1,
p*=1/γM²). MPI runs (`build-mpi/ns_solver` under `mpirun`) reproduce the
serial results bitwise or to round-off for every case below.

### 1. Oblique shock reflection (Euler) — `ObliqueShockRefEuler.toml`

Mach 2.9 flow over a flat wall; an oblique shock (β≈30°) enters through the
upper-left corner via prescribed inflow/post-shock states and reflects off
the bottom wall. Domain [0,3]×[0,1], 240×80 grid, NKFDS-MOVERS 1st order,
CFL 0.4.

States: inflow `(ρ,u,v,p)=(1, 2.9, 0, 1/1.4)`, post-shock
`(1.69997, 2.61934, −0.50633, 1.52819)` — both papers' §oblique-shock case.
Result: incident + reflected shock captured crisply at the theoretical
reflection point; wall pressure plateaus match oblique-shock theory
(p_ref ≈ 2 × p_post ≈ 3.06). Plot: `output/oblique/contours_p.png`.

### 2. Compression ramp in a channel (Euler) — `RampEuler.toml`

Mach 2 wind tunnel with a 15° ramp on the lower wall; leading-edge oblique
shock reflects from the top wall and interacts with the expansion fan at the
ramp shoulder. 240×80 grid, NKFDS-MOVERS-2O 2nd order, CFL 0.2.
Result: shock/reflection/expansion-fan triple structure; peak pressure 3.58
on the reflected-shock impingement. Plot: `output/ramp/contours_p.png`.

### 3. Forward-facing step (Euler, unsteady) — `FfsEuler.toml`

Woodward–Colella Mach 3 wind tunnel with a step, integrated to t = 4.0
(reached at iteration 3171 with CFL 0.6, NKFDS-MOVERS 1st order).
Result: lambda shock with triple point and slip stream, reflected shock /
expansion interaction — density plot `output/ffs/contours_rho.png`
(peak p ≈ 38 near the step corner).

### 4. Grid-aligned slip flow (Euler) — `SlipFlowEuler.toml` *(new)*

Mach 2 stream slipping under a Mach 3 stream at equal ρ and p on a 20×20
grid; tests exact capture of grid-aligned contact discontinuities.
Lower half of the left boundary prescribes `(1, 2, 0, 0.714286)`, upper half
`(1, 3, 0, 0.714286)`; all other boundaries transmissive.
Result: slip line held exactly — u jumps 2 → 3 across one cell with pressure
uniform to 10⁻³ everywhere after 5000 iterations. Plot:
`output/slip_flow/mach.png`.

### 5. Hypersonic flow past a half cylinder (Euler) — `HypersonicFlowEuler.toml`

Quirk-type carbuncle test: M∞ = 20 (and 6) flow onto a half-cylinder,
80×160 grid, NKFDS-MOVERS 1st order. The bow shock stands off cleanly with
no odd-even decoupling or carbuncle asymmetry. **Use a low CFL** — the
published CFL 0.4 goes unstable late in the run on this stiff case; CFL 0.02
(M=20) / 0.03 (M=6) run cleanly to 25k iterations with stagnation pressures
359 / 32.9 respectively. Plots: `output/hypersonic_m{20,6b}/mach.png`.

### 6. Blasius flat-plate boundary layer (NS) — `Blasius.toml`

Laminar M=0.5, Re=5000 flat plate; 160×64 stretched grid; ROE-2O +
Green-Gauss viscous fluxes. Validated against Blasius similarity theory:
δ99(x=1) = 0.073 vs analytical 0.0707, Cf within ~10% at x=1, smooth
monotone profile (see `output/swbli`-style analysis in git history).
Run: `... Blasius.toml BlasiusGrid16064.dat output/blasius 100000`.

### 7. Shock-wave/boundary-layer interaction (NS) — `SWBLI_NS.toml`

Hakkinen-style case: M=2.15, Re=10⁵ laminar flat plate with an impinging
oblique shock (β≈30.8°) prescribed through post-shock inflow states on
imin-top/jmax; 141×121 stretched grid; NKFDS-MOVERS-2O. Legacy-length run =
250k iterations (~30 min serial, faster with MPI).
Result: separation bubble x ∈ [0.78, 1.19] (height < 0.01), two-plateau wall
pressure 0.185 → 0.236, incident+reflected shock system. Plots in
`output/swbli/`; full discussion in `docs/Testcases.md`.

### 8. Supersonic viscous flow over a circular-arc bump (NS) — `BumpNS.toml` *(new)*

Channel [−1,2]×[0,1] with a 4%-thick arc bump (y = 0.04 sin πx, chord 1)
on the lower wall; M∞ = 1.4, Re = 8000. Grid generated by
`input/BumpGrid24080.dat` (240×80, 4.5% geometric stretching normal to the
wall; generator script inline in git history). Bottom: symmetry upstream of
the bump, no-slip adiabatic wall downstream; top: slip wall; supersonic
inlet/outlet. ROE-2O, CFL 0.25, 60k iterations.
Result: shock off the bump crest (p jumps 0.44 → 0.51 above it), supersonic
overshoot M=1.43, and a separated region x ∈ [0.79, 1.33] where the
reflected shock meets the boundary layer — the physics shown in Fig. 23–24
of Shrinath et al. Plots: `output/bump_ns/{contours_M,wall_pressure}.png`.

### Not yet configured (from the papers)

- **Double Mach reflection** (Shrinath §6.3.6): needs an internal initial
  condition (normal shock M=5.5 at x=0.25 inside the domain); the solver
  currently initialises uniform freestream only. Add restart-file support
  or an IC hook to enable.
- **Shock diffraction** (both papers): same IC requirement (moving shock
  initialised inside the domain), plus an L-shaped geometry that needs the
  cut-cell/cut-wall path exercised.

### Parity summary (this suite)

| Case | np=2 vs serial |
|---|---|
| Oblique shock / Ramp / FFS | inviscid j-slab halo — exact |
| Slip flow | exact (steady contact preserved) |
| Hypersonic half-cylinder | not MPI-tested (stiff; serial reference) |
| Blasius NS | ~2.5e-7 relative (round-off) |
| SWBLI NS | bitwise exact |
| Bump NS | same code path as SWBLI (gated BCs + interface gradients) |

