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
