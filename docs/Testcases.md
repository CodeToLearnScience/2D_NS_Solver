# Test Cases & User Guide

## Quick Start

```bash
# Build
cmake --preset release
cmake --build --preset release -j

# Run all 51 unit/parity tests
ctest --preset release --output-on-failure

# MPI build + tests
cmake -S . -B build-mpi -D CMAKE_BUILD_TYPE=Release -D NS_ENABLE_MPI=ON
cmake --build build-mpi -j
ctest --test-dir build-mpi -R mpi_halo --output-on-failure
```

## Running a Case

```bash
# Serial
./build/release/ns_solver configs/RampEuler.toml input/RampGrid24080.dat output/ 200

# MPI (4 ranks)
mpirun -np 4 build-mpi/ns_solver configs/RampEuler.toml input/RampGrid24080.dat output/ 200

# Restart from checkpoint
./build/release/ns_solver configs/RampEuler.toml input/RampGrid24080.dat output/ 100 \
    --restart output/checkpoint.rst
```

## Available Cases (configs/)

| Config | Grid | Physics | Scheme | BCs |
|---|---|---|---|---|
| `RampEuler.toml` | `RampGrid24080.dat` (240×80) | Euler, ND | NKFDS-MOVERS-2O | inflow, slip×2, transmissive |
| `Blasius.toml` | `BlasiusGrid16064.dat` (160×64) | NS, ND | ROE-2O | farfield, symmetry, no-slip, transmissive |
| `HypersonicFlowEuler.toml` | `HalfCylinderGrid80160.dat` (80×160) | Euler, ND | NKFDS-MOVERS | prescribed-inflow, transmissive×2, slip-wall |
| `FfsEuler.toml` | `FFSGrid24080.dat` (240×80) | Euler, ND | scheme-specific | per topology |
| `ObliqueShockRefEuler.toml` | `RampGrid24080.dat` | Euler, ND | KFDS-1O | per topology |
| `SWBLI_NS.toml` | `SWBLIGrid10565_2.dat` (105×65) | NS, ND | NKFDS-MOVERS-2O | per topology |

## Output Files

After running, the output directory contains:

| File | Format | Contents |
|---|---|---|
| `Residue_<case><nx>_<ny>.dat` | Tecplot ASCII | Iteration, Δt, time, L2(ρ), argmax, Cl, Cd, Cm |
| `<case>_IsoCont<nx>_<ny>Iter<N>.dat` | Tecplot ASCII | Node coordinates + rho, u, v, p, T, Mach |
| `<case>_IsoCont<nx>_<ny>Iter<N>.vtk` | VTK ASCII | StructuredGrid with density, pressure, velocity |

## Plotting Contours (ParaView)

1. Open the `.vtk` file in ParaView
2. Click **Apply** to load the dataset
3. Select **density**, **pressure**, or **velocity** from the coloring dropdown
4. For better visualization:
   - Toggle **Surface** → **Outline** or **Slice** for 2D data
   - Use **Calculator** filter for Mach number: `sqrt(vx^2+vy^2)/a`
   - Use **Contour** filter for iso-contour lines

Alternative — plot the `.dat` file with Python:

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Ramp_KFDS_maxU_24080_2nd_IsoCont240_80Iter200.dat", skiprows=3)
x, y = data[:,0], data[:,1]
rho, u, v, p, T, Mach = data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7]

nx, ny = 241, 81
xi = x.reshape(ny, nx); yi = y.reshape(ny, nx)

plt.figure(figsize=(12,4))
plt.subplot(131); plt.contourf(xi, yi, rho.reshape(ny,nx), 30, cmap='jet'); plt.colorbar(); plt.title('Density')
plt.subplot(132); plt.contourf(xi, yi, p.reshape(ny,nx), 30, cmap='jet'); plt.colorbar(); plt.title('Pressure')
plt.subplot(133); plt.contourf(xi, yi, Mach.reshape(ny,nx), 30, cmap='jet'); plt.colorbar(); plt.title('Mach')
plt.tight_layout(); plt.savefig("contours.png", dpi=150)
plt.show()
```

## Plotting Residual Convergence

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Residue_Ramp_KFDS_maxU_24080_2nd240_80.dat", skiprows=2)
iters, dt, t, res = data[:,0], data[:,1], data[:,2], data[:,3]

plt.semilogy(iters, res)
plt.xlabel("Iteration"); plt.ylabel(r"$\|\Delta\rho\|_2$")
plt.title("Residual Convergence")
plt.grid(True)
plt.savefig("residual.png", dpi=150)
plt.show()
```

## Unit / Parity Test Suite (51 tests)

Run: `ctest --preset release --output-on-failure`

### Config Layer (`test_config.cpp`) — 16 tests

| # | Test | What it validates |
|---|---|---|
| 1 | MinimalValidConfigParsesWithDefaults | TOML parse; legacy defaults (γ=1.4, Pr=0.72, R=286.92) applied |
| 2 | MissingFileIsCleanError | Clean error for nonexistent config file |
| 3 | MissingRequiredKeyIsReported | Missing `numerics.cfl` caught with key path |
| 4 | UnknownKeyIsRejected | Typo'd keys rejected (not silently ignored) |
| 5 | InvalidEnumValueListsAllowedOptions | `"warp-drive"` rejected; valid options listed in error |
| 6 | MalformedTomlReportsParseError | Syntax error reported with file+line |
| 7–15 | ConfigValidate.* | Cross-field validation: steady needs max_iter, unsteady needs total_time, restart needs file, ND needs reference_state, isothermal needs T_wall, segment bounds checked, multiple errors batched |
| 16 | ConfigRoundTrip.ToTomlSurvivesReload | Serialize→parse round-trip preserves all fields |

### Fields & Ghost Model (`test_field.cpp`) — 5 tests

| # | Test | What it validates |
|---|---|---|
| 17 | GhostAwareIndexing | Negative-index ghost access; at() throws OOB; assert-guarded operator() |
| 18 | ContiguousJFastestMatchingLegacy | Memory layout matches legacy double** (j-fastest rows) |
| 19 | FillAndGhostsCopy | fill_ghosts_copy() reproduces legacy convention exactly |
| 20 | MdspanViewCoversWholeCanvas | View2D[i,j] covers full ghosted canvas |
| 21 | PlanesAreContiguousAndIndependent | MultiField planes are contiguous, independent |

### Mesh & Metrics (`test_mesh.cpp`) — 6 tests

| # | Test | What it validates |
|---|---|---|
| 22 | ReadsLegacyIxFastestFormat | Legacy grid reader (i-fastest, comment-skip) |
| 23 | AppliesScaling | grid.scaling factor multiplies coordinates |
| 24 | RejectsTruncatedAndOverlongFiles | Clean errors for malformed grids |
| 25 | UniformCartesianGrid | Face vectors (∓1,0)/(0,∓1), area=1 on Cartesian |
| 26 | AffineCellAreaEqualsJacobianDet | area = \|det J\| for affine map |
| 27 | GeneralQuadMatchesShoelace | Skewed quad vs independent shoelace formula |
| 28 | HalfCylinder45LoadsAndHasPositiveArea | Real grid: all areas > 0, exact node coords |

### IO (`test_io.cpp`) — 3 tests

| # | Test | What it validates |
|---|---|---|
| 29 | EmitsLegacyCompatibleStructuredGrid | VTK header, node ordering, POINT_DATA blocks |
| 30 | RoundTripIsBitExact | Binary restart write→read via memcmp |
| 31–32 | RejectsNonRestartFiles / TruncatedPayload | Corruption detection |

### Physics Parity (`test_physics.cpp`) — 6 tests

| # | Test | What it validates |
|---|---|---|
| 33 | GreenGaussFaceGradientsMatchLegacyBitwise | 42 golden values from legacy objects — exact match |
| 34 | ViscousFluxMatchesLegacyBitwise | 9 golden values — exact |
| 35 | DependentVariablesMatchLegacyBothFormulations | 48 values across ND and DIM — exact |
| 36 | PrimConsRoundTripIsIdentity | EOS prim→cons→prim identity |
| 37 | AnalyticRelationsHold | Sound speed, temperature, pressure from E |
| 38 | DimensionalFormEqualsClassicSutherland | μ(T) = C₂T^{1.5}/(T+S) recovered |

### Numerics Parity (`test_numerics.cpp`) — 10 tests

| # | Test | What it validates |
|---|---|---|
| 39 | LlfFirstOrderMatchesLegacyBitwise | Complete LLF chain bitwise vs legacy goldens |
| 40 | LlfSecondOrderMatchesLegacyBitwise | MUSCL LLF chain bitwise |
| 41 | NkfdsMoversFirstOrderScheme7MatchesLegacyBitwise | NKFDS-MOVERS 1st order chain bitwise |
| 42 | RoeSecondOrderScheme24MatchesLegacyBitwise | Two-wave Roe + VanAlbada MUSCL chain bitwise |
| 43 | NkfdsMoversSecondOrderScheme30MatchesLegacyBitwise | NKFDS-MOVERS 2nd order chain bitwise |
| 44 | ConstantFieldGivesZeroRhsBothOrders | Free-stream preservation (diss exactly zero) |
| 45 | MinmodFluxLimiterBoundsAndValues | Limiter clamps [0,1] and known values |

### Decomposition (`test_parallel.cpp`) — 3 tests

| # | Test | What it validates |
|---|---|---|
| 46 | RangesPartitionGlobalRows | Contiguous, non-overlapping slabs cover global domain |
| 47 | OwnershipFlags | jmin/jmax/lower/upper neighbor flags correct |
| 48 | RejectsInvalid | More ranks than rows; out-of-range rank |

### Restart Continuation (`test_restart.cpp`) — 1 test

| # | Test | What it validates |
|---|---|---|
| 49 | SerialBitwiseContinuation | Checkpoint→load→continue produces identical state to uninterrupted run |

### Physics Verification (`test_physics_verification.cpp`) — 2 tests

| # | Test | What it validates |
|---|---|---|
| 50 | FreeStreamPreservationAcrossSchemes | Uniform field → exactly zero dissipation (LLF, ROE-TV, H1, H2) |
| 51 | PositivePressureAfterShockTubeIteration | 200 iterations of ramp-Euler produce no NaN/negative states |

## MPI Tests (`test_mpi_halo.cpp`) — 2 tests

Run under `mpirun`: `ctest --test-dir build-mpi -R mpi_halo`

| Test | Ranks | What it validates |
|---|---|---|
| mpi_halo_np2 | 2 | Halo exchange: ghost cells match neighbour real rows exactly |
| mpi_halo_np3 | 3 | Same with interior rank (two interfaces) |

## Known Limitations

| Issue | Impact | Status |
|---|---|---|
| Parallel HDF5 IO not implemented | Restart is root-only binary format; no parallel IO | Deferred |
| METIS 2D decomposition not implemented | Only 1-D j-slab decomposition available | Low priority for current use |
| CI/CD pipeline not configured | Tests run manually only | Enhancement |
| Surface pressure distribution writer not ported | Cl/Cd/Cm in residue file work; surface distribution .dat missing | Low priority |

**Fixed**: NS viscous chain crash was caused by `green_gauss_face_gradients` accumulating
into `gradfi`/`gradfj` without zeroing them first. Gradients grew unboundedly each iteration.
Fix: call `.fill(0.0)` before accumulation (matching legacy `Initialize_Gradients`).

## Benchmarking

```bash
scripts/bench.sh 300
```

Example output (ramp-Euler, 300 iterations):

```
legacy                         2.69 s  ( 111.6 iter/s)
new serial                     2.15 s  (  139.4 iter/s)   ← ~25% faster
new mpi -np 4                  1.03 s  (  292.1 iter/s)   ← ~2.6× faster
```
