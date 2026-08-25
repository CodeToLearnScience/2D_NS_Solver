# Refactoring Plan: 2D NS Solver → Modern C++ / MPI / CMake

Status: **living document**. Update phase tables as work lands.

---

## 1. Goals

1. Rewrite the solver core in modern C++ (C++20, C++23 where toolchain allows).
2. Replace the positional `.cfg` + `GridTop*.dat` inputs with a validated **TOML** config
   (JSON optional backend).
3. **MPI-parallelize** via domain decomposition of the structured grid.
4. Migrate the build from the hand-rolled Makefile to **CMake** (presets, FetchContent,
   warnings-as-errors, sanitizers, CI-ready).
5. Establish **Google Test** based unit + verification + regression testing.

## 2. Current state (baseline audit)

~14,400 LOC, 47 `.cpp` files in `src/`, headers in `inc/`, built with `g++ -Wall -O3` Makefile.
Flow: `main.cpp` → positional input parse → read grid + separate topology file → allocate ~30
global arrays → march `Solver()` (SSPRK1–3, local Δt) calling scheme-specific `Diss_*`
inviscid-flux routines, Green-Gauss viscous fluxes for NS, segment-table boundary conditions,
Tecplot/VTK output, restart support.

Schemes present: LLF, Roe (1st/2nd order), MOVERS variants (H1/H2/LE1/LE2/NWSC), KFDS,
NKFDS-MOVERS, ECCS. BC flags: 10 prescribed inflow, 20 transmissive, 30 Euler wall,
40 NS adiabatic wall, 50 NS isothermal wall, 60 far-field, 70 symmetry, 80 cut wall.

### 2.1 Structural problems

| Problem | Where |
|---|---|
| ~60 mutable globals threaded through every function | `inc/global_declarations.h`, `src/initialize_variables.cpp` |
| Raw `double***` jagged arrays; manual new/delete; leaks | everywhere; e.g. `main.cpp` frees one pointer array via comma operator |
| Variable-major layout `cv[var][i][j]` with j-outer/i-inner loops → strided access | all hot kernels |
| Cryptic index zoo (`id1/id2/jd1/jd2/imax`) for ghost cells | every file |
| Scheme/BC dispatch by bare int codes in if-chains | `solver.cpp`, `boundary_conditions.cpp` |
| Hard-coded relative paths `../input/`, `../output/` | `read_grid.cpp`, `read_grid_top.cpp`, `main.cpp` |
| Fragile positional parsing by repeated `getline`+`>>` | `read_solver_input.cpp`, `read_grid_top.cpp` |

### 2.2 Bugs found and fixed (2026-08, pre-baseline)

These were corrected **before** capturing regression baselines so that goldens reflect fixed
behavior. Any legacy results produced earlier are not comparable byte-for-byte.

| Bug | Location | Fix |
|---|---|---|
| Face-normal y-component computed as `sx/ds` instead of `sy/ds` → wrong normal velocity in spectral radii → wrong Δt | `src/time_step.cpp`, both `Time_Step_Euler` and `Time_Step_NS` | `ny = sy / ds;` (4 sites). Verified against metric convention in `grid_computations.cpp`: `si[0]=-Δy` (x-comp), `si[1]=+Δx` (y-comp) |
| `MinModLim` missing return path when `|Ur|==|Ul|` (UB) | `inc/basic_functions.h` | Rewrote as guard-on-opposite-signs then smaller-magnitude selection |
| No-op string comparison `time_step == "global";` | `src/time_step.cpp` (both functions) | Removed |
| Uninitialized `cs` read in error-diagnostic print | `src/time_step.cpp` | Initialized to 0.0 |
| `while(!infile.eof())` wrapping entire grid read loop → re-read pass zeroes first grid point when file ends with newline | `src/read_grid.cpp` | Single guarded pass with failure diagnostics |
| Scheme banner printed from `time_accuracy` instead of `inviscid_scheme` | `src/read_solver_input.cpp` | Corrected variable |
| `Log_Mean` else-branch inverted `f/ln(ξ)` and lost the sign → **negative** logarithmic mean for density jumps beyond ~±10% (e.g. −1.56 instead of +1.4427 for states (2,1)) | `inc/basic_functions.h` | Broken branch replaced with exact `(al−ar)/ln(al/ar)`; small-jump series branch untouched. NOTE: currently zero callers in `src/` — becomes live when Roe/KFDS paths are ported |
| `get_str_between_two_str` searched the closing delimiter from string start, finding the opening one again → returned everything *after* delimiter #1 instead of the text *between* | `inc/basic_functions.h` | Closing search starts past the opening delimiter; missing-closer falls back to end-of-string. Existing restart-file callers unaffected (they parse leading numerics) |

Verification findings from baseline capture:

- The committed 80×160 hypersonic half-cylinder case runs cleanly post-fix
  (200-iteration smoke: residuals finite, monotonically decreasing).
- The 45×45 half-cylinder grid is **not viable at Mach 20** with these input parameters:
  pristine pre-fix code survives only 3 iterations (masked by the wrong Δt), fixed code hits
  the negative-pressure guard at iteration 1. Treated as a case-definition issue, not a code
  regression; excluded from baselines.

### 2.3 Legacy input-file inconsistencies (found during Phase 2)

The positional reader in `read_solver_input.cpp` silently mis-parses files whose field layout
drifted from what it expects (`>> int` on `2.15` yields `2`; missing blocks shift every later
field). Four committed inputs are affected and are **not** convertible:

| File | Problem |
|---|---|
| `InputHypersonicFlowNS.dat` | `space_accuracy = 2`, which the legacy solver itself rejects at runtime |
| `InputSWBLI_NS_Dimen.dat` | extra/misordered numeric fields; reader would consume garbage values |
| `InputSWBLI_NS_ND1.dat` | old format predating the `disp_freq`/`outp_freq`/`tot_time` blocks (commit `91da8d7`) |
| `InputWedgeReflection.dat` | same old format; e.g. `outp_freq` would read `257354` |

If any of these cases is still needed, recreate its input from a working template — the new
TOML loader would have rejected them loudly instead of running with corrupted settings.

Not-yet-reviewed: numerics inside `diss_roe*.cpp`, `diss_movers*`, `ECCS.cpp`, `bc_farfiled.cpp`
(≈5,500 LOC). Deep review is deferred to Phase 5 where each scheme is ported against golden
outputs rather than eyeballed.

Additional legacy issues found while porting the physics (Phase 4), fixed by the rewrite:

- `Flux_Viscous` leaks its `fv = new double[3]` scratch buffer on **every call**.
- `Gradient_FaceI/J` close with scalar `delete fgx;` on `new double[3]` memory
  (alloc/dealloc mismatch, UB; ASan flags it).
- `Flux_Viscous` writes through the global `diss` pointer — the class of global-state
  coupling this refactor removes.

Parity methodology: a probe program is compiled against the legacy objects and prints
`%.17g` outputs for `Dependent_Variables_One`, `Gradient_FaceI/J`, and `Flux_Viscous` on
the RampGrid6020 geometry with a deterministic analytic field set. The values are checked in
as `tests/goldens_physics.inc`; the new kernels must reproduce them exactly
(`EXPECT_DOUBLE_EQ`). The transcription preserves legacy arithmetic order per element,
including quirks such as the un-normalized top row of `gradfi` / right column of `gradfj`
(never read by consumers) and the differing quarter-average orderings between the I and J
gradient passes.

## 3. Target architecture

```
2D_NS_Solver/
├── CMakeLists.txt              # root; ns_core lib + apps + tests
├── CMakePresets.json           # gcc/clang × debug/release/asan
├── REFACTORING_PLAN.md         # this document
├── configs/                    # migrated TOML cases (+ tools/cfg2toml output)
├── src/
│   ├── main.cpp                # thin driver: args → config → run
│   ├── config/                 # TOML loader → typed Config structs + validation
│   ├── mesh/                   # StructuredMesh: nodes, metrics (face vectors, areas)
│   ├── fields/                 # Field<T>: contiguous vector + ghost-aware View2D
│   ├── physics/                # IdealGas EOS, prim↔cons, Sutherland viscosity, gradients
│   ├── numerics/               # inviscid fluxes (LLF/Roe/MOVERS/KFDS/ECCS), limiters, MUSCL
│   ├── bc/                     # enum class BcType + per-type appliers on segments
│   ├── time/                   # Steppers: ForwardEuler, SSPRK2, SSPRK3
│   ├── parallel/               # decomposition, halo exchange, reductions (MPI)
│   ├── io/                     # grid readers, restart (binary), VTK/XDMF, surface writers
│   └── solver/                 # driver loop orchestration
├── tests/                      # Google Test unit + verification + MPI regression
├── scripts/                    # capture_regression.sh, compare_regression.sh
├── tools/cfg2toml/             # legacy config migration tool
└── legacy/                     # (eventually) frozen pre-refactor sources, kept until parity
```

### 3.1 Key design decisions

- **No globals**: a `SimulationContext` owns config/mesh/fields; kernels take views.
- **Fields**: contiguous `std::vector<T>` with explicit ghost layers (`ng`, default 3),
  exposed through `ns::View2D<T>` — a minimal 2-D view whose `[i, j]` subscript and
  `extent()` API deliberately mirror `std::mdspan`. Deviation from the original plan:
  libstdc++ (checked on GCC 15) still lacks `<mdspan>`, so the real type is a drop-in
  swap once available. Bounds-checked `at()` accessors complement assert-guarded
  `operator()` (zero cost in release hot loops). `fill_ghosts_copy()` reproduces the
  legacy ghost convention exactly (each ring copies the nearest interior layer).
- **Concepts**: `template<Equation E, Flux F>` constraints; runtime selection via enum +
  factory at step granularity (dispatch cost negligible vs stencil sweeps).
- **BC framework**: `enum class BcType`, `BoundarySegment {type, edge, start, end, state}`;
  topology moves *into* the TOML file.
- **Fallible ops**: `std::expected<T, Error>` for parse/load paths — no `exit()` mid-parse.
- **C++23 features used**: `std::expected`, `std::format`, multidim-subscript views,
  deducing-this where helpful. C++20 baseline otherwise (`ranges`, concepts, spaceships).

## 4. Configuration design (TOML)

Library: **toml++** (header-only, MIT). Optional JSON mirror via nlohmann/json later if needed.

```toml
[case]
name = "HyperCylM20"
equations = "euler"             # euler | navier-stokes
formulation = "nondimensional"  # nondimensional | dimensional
flow = "steady"                 # steady | unsteady

[grid]
file = "HalfCylinderGrid80160.dat"
scaling = 0.0381

[numerics]
inviscid_scheme = "nkfds-movers" # llf | roe | movers-h1 | movers-le1 | kfds | eccs | ...
order = 1                        # 1 | 2
limiter = "minmod"               # minmod | van-albada | venkatakrishnan
time_method = "ssprk3"           # forward-euler | ssprk2 | ssprk3
cfl = 0.4
residual_smoothing_eps = 0.0     # epsirs

[physics]
gamma = 1.4                      # hard-coded today; becomes configurable
Re_inf = 257354.0
Mach_inf = 20.0
alpha_deg = 0.0
p_inf = 101325.0
T_inf = 298.0

[physics.reference_state]        # nondimensional runs only
pressure = 100000
temperature = 300
density = 1.12
velocity = 1

[run]
max_iterations = 100000
total_time = 30.0                # unsteady stop time
restart = false
restart_file = ""

[output]
directory = "output"
display_frequency = 1000
write_frequency = 10000
formats = ["vtk", "surface", "residual"]

[[boundary]]
edge = "jmin"                    # imin | imax | jmin | jmax
start = 2                        # node range along the edge
end   = 161
type  = "farfield"

[[boundary]]
edge = "imax"
start = 2
end   = 161
type  = "euler-wall"
```

Loader validates everything up-front with `file:line` diagnostics (toml++ node sources) and a
schema check (unknown keys rejected). `tools/cfg2toml` converts each legacy pair
(`Input*.cfg` + `GridTop*.dat`) into this format so all existing cases migrate mechanically.

## 5. MPI design

| Aspect | Decision |
|---|---|
| Decomposition | 1-D slabs along `j` first (N−1 interfaces); upgradeable to 2-D bisection/METIS |
| Ghost depth | `ng = 3` (covers 2nd-order MUSCL + Green-Gauss face gradients) |
| Halo exchange | Pack contiguous per-variable planes → `MPI_Sendrecv`; derived datatypes only after profiling |
| Ordering | Halo fill **before** BC application each stage/substage (encode in driver) |
| Reductions | Residual norms, min Δt, forces via `MPI_Allreduce` |
| I/O | Phase A: rank-local writes + aggregator; Phase B option: parallel HDF5 + XDMF |
| Restart | Fixed-width binary, typed header, MPI-agnostic layout |
| Reproducibility | Reduction order ⇒ results depend on rank count within FP tolerance; documented + tested |

## 6. Build system (Phase 1 — done)

- Targets: `ns_legacy` (frozen original code, output name `2D_NS_Solver`), later `ns_core` +
  `ns_solver`; `ns_tests` (GTest); `cfg2toml`.
- Options: `NS_BUILD_TESTS` (ON), `NS_ENABLE_WARNINGS_AS_ERRORS` (OFF), sanitizers preset.
- GTest: `find_package` first, FetchContent fallback (pinned release).
- Presets in `CMakePresets.json`; ccache-friendly; Ninja or Makefiles both fine.
- Old `Makefile` retained until parity is proven, then removed.

## 7. Testing strategy (Google Test)

| Layer | Content |
|---|---|
| Unit | EOS round-trips, sound speed; flux constancy on uniform states; Roe eigenvalues; limiter properties (min-max, symmetry); minmod incl. tie case (regression for the UB fix); quad metrics on known geometry; config parse/validation errors; halo pack/unpack logic |
| Verification | Sod problem (1D-in-2D); isentropic vortex; order-of-accuracy convergence (L1 rates ≈1 and ≈2); laminar flat plate vs `input/Blasius*.dat` references |
| Regression | Golden outputs vs legacy solver for ramp-Euler, Blasius-NS, hypersonic half-cylinder cases captured by `scripts/capture_regression.sh` (byte-compare via SHA256 with `scripts/compare_regression.sh`) |
| MPI | Same case on 1/2/4 ranks; residual trajectories match within tolerance; smoke scaling test under `mpirun` |

## 8. Phased roadmap

| Phase | Scope | Gate |
|---|---|---|
| 0 ✅ | Baseline audit + bug fixes (see §2.2) | fixes reviewed; 11/11 unit tests green |
| 1 ✅ | CMake + GTest scaffolding + regression scripts | `ctest` green; goldens captured post-fix for ramp-Euler, Blasius-NS, hypercyl-80160 (200-iter smokes); compare tool idempotent |
| 2 ✅ | TOML config loader + validator + cfg2toml tool; unit tests | 33/33 tests green; all 6 convertible legacy cases migrated to `configs/*.toml` via `scripts/migrate_configs.sh` and re-validated (4 known-incompatible inputs documented in §2.3) |
| 3 ✅ | Mesh/fields/io infra (ghost-aware fields, metrics, legacy grid reader, VTK writer, binary restart) | 49 tests green incl. hand-computed geometry (uniform/affine/skewed-quad shoelace), real HalfCylinder grid (all areas > 0), VTK content, restart bit-exact round-trip, truncation/corruption rejection; ASan/UBSan clean |
| 4 ✅ | Physics port: EOS, Sutherland transport, dependent variables, Green-Gauss face gradients, viscous fluxes | **bit-for-bit parity** with legacy kernels via golden values generated from the legacy objects (gradients: 42 values, viscous flux: 9, dependent vars: 48 across both formulations); plus EOS round-trip / classic-Sutherland property tests. 55 tests green; ASan/UBSan clean |
| 5 ✅ | Inviscid scheme ports. **5a:** shared face-flux infrastructure, LLF 1st/2nd-order. **5b:** NKFDS-MOVERS 1st-order (scheme 7), two-wave Roe w/ Van Albada MUSCL (24), NKFDS-MOVERS 2nd-order (30) — the three schemes used by the regression baselines. All verified bitwise against complete legacy chains; legacy quirks preserved and documented (third `MaxEigVal` variant without `+a`; `Movers()` rounds jumps to 3 decimals; Mahalanobis entropy switch driven by global `Cp = 1/((γ−1)M∞²)` in ND runs). Remaining legacy-only variants (MOVERS H1/H2/LE1/LE2/NWSC, KFDS/ECCS, Roe-TV-1st) ported on demand — none are used by any migrated config | 62 tests green incl. 5 scheme-parity chains; ASan/UBSan clean |
| 6 ✅ | Steppers, BC framework, driver → single-rank end-to-end (`ns_solver`, ramp-Euler case: scheme 30 + forward Euler + prescribed-inflow/slip-wall/transmissive) | Residual trajectory ratio 1.000006 @ 200 iterations; field agreement ~1e-4 relative (chaos floor, see below); first-iteration state byte-aligned; 62 unit/parity tests green; ASan/UBSan clean |

**Phase-6 end-to-end parity notes.** Bit-for-bit file equality holds for the entire
first iteration; from iteration ~3 onward the legacy scheme's *quantized* wave-speed
estimator (`Movers()` rounds flux/state jumps to 3 decimals) turns sub-roundoff
codegen-order differences into discrete branch flips, which chaotically amplify to
~1e-4 relative field differences by iteration 200 while the convergence trajectory
stays identical (final-residual ratio 1.000006). Kernel-level bitwise gates remain
the Phase-4/5 goldens; cross-codegen bitwise identity of the full driver is not an
achievable target for this discretization. Reproduce via
`scripts/check_e2e_parity.sh [iters]`.

Additional legacy quirks discovered & replicated in Phase 6 (see §2.2 for the rest):
- Startup order matters: legacy runs `Dependent_Variables` **before** the pre-loop
  boundary pass, so all planes start from cv-reconstructed values, one ULP off the
  analytic freestream.
- `Mahalanobis` reads a **separate lowercase global `cp`** that is zero-initialized,
  never set by the initialization, and clobbered by `Forces()` with a local pressure
  coefficient after every residue evaluation -- the entropy switch therefore runs on
  iteratively-mutated garbage. The port carries this as an explicit `mah_cp` state.
- Legacy zeroes `diss` only over interior cells each step (ghost-ring dissipation
  accumulates stale scatter); harmless because the RHS seeding only reads interiors,
  but replicated by construction.
| 7 ✅ | MPI: 1-D j-slab decomposition, packed `MPI_Sendrecv` halo exchange, global min-Δt/residual/force reductions, root-aggregated writers, gather-based solution output. Gate (ramp-Euler @200 iters): Δt bitwise-identical across 1/2/4 ranks (order-free MIN); final residual matches legacy to ~1e-5 abs for np=1/2 and ~2e-4 for np=4; fields within ~0.5% column-relative across rank counts — dominated by the documented `mah_cp` quirk whose "last integrated face" is decomposition-dependent (inherent to replicating the legacy scheme; see §6 notes). Distributed binary restart deferred to Phase 8 alongside parallel-HDF5 I/O. Tests: decomposition partitioning (serial), halo exchange correctness (`mpirun -np 2/3`), full e2e at three rank counts via `scripts/check_mpi_parity.sh`; serial suite still 65/65 green |
| 8 ✅ | Perf/polish: benchmark harness (`scripts/bench.sh`: new serial ~25% faster than legacy on ramp-Euler; np=4 ≈2.6× legacy on this small case), restart continuation (write/load, serial bitwise-continuation test, `--restart` flag incl. MPI), `NS_ENABLE_NATIVE_ARCH` option (documented FMA/parity caveat), ASan/UBSan-clean suites, README rewrite. Deferred: parallel-HDF5 IO, METIS 2-D decomposition, clang-tidy (not installed) | benchmarks recorded; 67 tests green; e2e np=1/2/4 consistent |

## 8.1 Completion & cutover plan (Phases 9–11)

User decision: **port everything** (option a); fill legacy gaps with correct implementations;
physical removal of the legacy tree after validation (recoverable via git history). Modern-
ization freedom granted — physics correctness is the bar, not exact transcription.

| Phase | Scope | Gate |
|---|---|---|
| 9 | **Coverage completion.** Remaining inviscid schemes ported under golden parity: MOVERS(1)/H1(2)/LE1(3) [Batch A], ROE-TV(5)+KFDS(6) [B], MOVERS-H2(22)/LE2(23)/ROE-TV-2O(25)/KFDS-2O(26) [C], ECCS(27)/MOVERS-NWSC(29) [D]. **Gap fills:** scheme-21 (2nd-order H2 combination) implemented properly; entropy-switch `cp` defect fixed (physically consistent heat-capacity model replaces the Forces()-clobbered global; documented trajectory shift). **Driver completion:** viscous chain wired for NS, dimensional formulation, no-slip adiabatic/isothermal walls + characteristic far-field + cut BCs, real SSPRK2/SSPRK3 steppers, unsteady stop-on-total_time, surface/forces writer | all six migrated configs e2e serial+MPI ≥100 iters; bitwise goldens per batch; Blasius-NS residue/surface vs frozen baselines |
| 10 | **Physics verification suite** (legacy-independent): Sod & Lax shock tubes (L1 errors, discontinuity positions), isentropic vortex, order verification (~1st/~2nd), flat-plate skin friction vs Cf≈0.664/√Reₓ, free-stream preservation across every registered scheme | physics suite green; convergence orders within tolerance |
| 11 | **Golden freeze → validation → deletion.** Regenerate & commit all goldens/baselines (incl. surface outputs) before removal; feature-parity sweep (every convertible input has TOML coverage); cutover checklist (configs e2e np=1/2/4 · suite green · sanitizers clean · trajectories within floors); delete legacy `.cpp` bodies, `inc/`, `Makefile`, `ns_legacy` target, Riemann driver, obsolete `input/*.cfg`+`GridTop*.dat` (grids kept); docs rewrite | checklist complete; working tree contains only modern stack |

Post-cutover policy: frozen bitwise goldens become tolerance-based regression pins; the
physics suite is the primary correctness net. Upstream-dead code (`Diss_MOVERS2_Prim`,
`Diss_ZB2`, scheme-21 stub, second conflicting `Diss_Roe1`) is not ported; config
validation rejects those enums with clear errors.

## 9. Dependencies

| Package | Version pin | Purpose |
|---|---|---|
| Google Test | v1.15.x (FetchContent) | unit/regression tests |
| toml++ | v3.4.0 (FetchContent) | config parsing |
| MPI | system (OpenMPI present on dev box) | parallelism |
| fmt / CLI11 | optional | logging / argv (decide at Phase 2) |

## 10. Open decisions

1. ~~TOML-only vs dual TOML/JSON~~ → default TOML-only; JSON backend only if a concrete need appears.
2. Compiler floor: GCC ≥ 13 recommended (`std::expected`/`std::format` completeness); swap `View2D` → `std::mdspan` when libstdc++ ships it.
3. Whether historical published results need bitwise reproduction (would require keeping the old
   Δt bug behind a compat flag) — presumed no, given fixes were requested.

## 11. How to build & test (current state)

```bash
cmake --preset release          # or: debug / asan
cmake --build --preset release -j
ctest --preset release          # unit tests
scripts/capture_regression.sh   # rebuilds, runs shortened legacy cases, archives goldens
scripts/compare_regression.sh <dirA> <dirB>   # SHA256 byte comparison
```

Note: cmake was installed locally at `~/.local/opt/cmake-3.31/bin` on the dev box
(no sudo available); add it to PATH or install distro-wide when possible.
