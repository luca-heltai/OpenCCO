# Local macOS build with Spack dependencies

Requires CMake 3.16 or newer, Ninja, a C++17 compiler, and Eigen, Boost and Zlib.
Ceres 2.0.0 is bundled; DGtal sources must be available in `ext/DGtal`.
The validated DGtal revision is `cc6e2132d8c26ba11b303bdf85efa3b3937942a5`
(1.5.beta). Do not overwrite an existing dependency checkout to change revisions.

The commands below run from the OpenCCO root. On the validated machine, Spack
is available at `/Applications/deal.II.app/Contents/Resources/spack/bin/spack`;
put that directory on `PATH` first if necessary. Inspect `spack find -p` before
selecting packages. This build used Eigen 3.4.1, Boost 1.83.0 and Zlib-ng 2.3.3
with its Zlib compatibility interface.

```sh
spack find -p
EIGEN=$(spack location -i eigen)
BOOST=$(spack location -i boost)
ZLIB=$(spack location -i zlib-ng)
PREFIX="$PWD/build/deps/install"

# Apply once to pristine dependency sources. For an already-patched tree,
# use patch --dry-run -R to verify the patch is present instead of reapplying.
patch -p1 -d ext/ceres-solver-2.0.0 < cmake/patches/ceres-cmake.patch
patch -p1 -d ext/DGtal < cmake/patches/dgtal-build.patch

cmake -S ext/ceres-solver-2.0.0 -B build/deps/ceres -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DEigen3_DIR="$EIGEN/share/eigen3/cmake" \
  -DMINIGLOG=ON -DMINIGLOG_MAX_LOG_LEVEL=0 -DGFLAGS=OFF -DSUITESPARSE=OFF -DCXSPARSE=OFF \
  -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DBUILD_BENCHMARKS=OFF
# Level0 suppresses verbose traces, retaining INFO/warnings/errors/fatal logs.
cmake --build build/deps/ceres --parallel 4
cmake --install build/deps/ceres

cmake -S ext/DGtal -B build/deps/dgtal -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DEIGEN_ROOT="$EIGEN/include/eigen3" \
  -DCMAKE_PREFIX_PATH="$BOOST;$ZLIB" \
  -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF -DBUILD_POLYSCOPE_EXAMPLES=OFF
cmake --build build/deps/dgtal --parallel 4
cmake --install build/deps/dgtal

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PREFIX;$BOOST;$ZLIB" \
  -DCeres_DIR="$PREFIX/lib/cmake/Ceres" \
  -DDGtal_DIR="$PREFIX/lib/DGtal" \
  -DEigen3_DIR="$EIGEN/share/eigen3/cmake" \
  -DCMAKE_CXX_COMPILER_LAUNCHER= -DCMAKE_C_COMPILER_LAUNCHER=
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure --no-tests=error
```

Use the same Eigen version for Ceres and OpenCCO; Ceres checks the exact version.
The dependency prefix is local to `build/`: the old `ext_inst/` installation and
Spack packages are not modified. Compiler launchers are disabled because the
sandbox cannot write to the default `~/.ccache` directory. A writable ccache
installation/cache can be enabled separately.

The dependency patches raise Ceres' obsolete CMake compatibility version, use
Boost's CMake package instead of deprecated FindBoost, and make DGtal's unused
format/ccache helper downloads opt-in. They also report missing optional Graphviz
as a status message rather than a build warning. Warning diagnostics are not
suppressed. Optional DGtal GUI, CGAL, GMP and HDF5 features are not enabled;
headless 2D/3D tree generation and mesh export do not require them.

## Validation

All executable targets compile. CTest registers the four self-contained legacy
tests plus native-input, domain and generation-CLI regressions. Legacy tests
have isolated output directories and timeouts. Printed `KO` and
`TEST ERRORS` markers count as failures, even when a test exits zero. These legacy
tests are not a comprehensive numerical verification of the physical model.
`testCompCCO` is still compiled but is not registered: it needs the external
`../Data/Nt10_kt10_s42_M301/NoCom_Nt10_s42_M301_` reference files
(`TerminalSeeds.txt`, `distal.txt`, `proximal.txt`), which are not distributed.
See fork issue #2; do not fabricate reference data or count this program as passed.

Run small generator checks in scratch directories to avoid cluttering the source:

```sh
ROOT="$PWD"
mkdir -p build/smoke/2d build/smoke/3d
(cd build/smoke/2d && "$ROOT/build/bin/generateTree2D" -n 20 -a 20000 -x tree.xml)
(cd build/smoke/3d && "$ROOT/build/bin/generateTree3D" -n 20 -a 20000 -x tree.xml)
```

For both implicit sphere/circle and square (`-s`) domains, the validated outputs
contain 40 XML nodes and 39 edges for 20 terminals, with finite numeric values.
SVG/EPS and OFF exports were also checked. A 2D mask run using
`-d "$ROOT/Samples/shape3.pgm"` passed the same structural checks.

The bundled 3D bunny mask has restrictive geometry. In the original build,
a 20-terminal, area-20000 Release run generated only eight terminals yet returned
zero; a Debug run exceeded a 90-second timeout. Current generators return
nonzero when their finite search budget leaves a partial tree, while retaining
partial exports and reporting the actual/requested counts. A seed-23 bounded
run at the README root retained one terminal and correctly returned 1.
These are different sampling policies, not a speed or trajectory comparison.
Partial output is not a successful full-tree test. See fork issues #3 and #12.

## Finite growth runs

Both generators accept `--seed`, `--searchTrials` (default 100 candidate searches
per insertion) and `--samplingTrials` (default 1000 total samples per search).
Counts and budgets must be positive. The seed defaults to wall-clock time;
set it explicitly in recorded experiments:

```sh
./build/bin/generateTree3D -n 8 --seed 23 --searchTrials 100 --samplingTrials 1000
```

Spacing starts at the native distance threshold and retains the legacy 0.9
relaxation factor after each 100 failed samples, now within a finite **total**
sampling budget. It no longer relaxes indefinitely. The earlier attempted
fixed-threshold-only policy prevented ordinary startup and was rejected.
Each sample invokes the domain controller; a custom controller must itself
terminate. Count budgets are not wall-time guarantees, so use an external
time/resource limit for long runs.

`tryGenerateNewLocation` reports budget exhaustion without returning a usable
candidate. The legacy point-returning `generateNewLocation` now throws
`std::runtime_error` on exhaustion. `expandTree` returns a completion boolean,
stops on an exhausted insertion budget, and counts success even on the last
allowed attempt. Callers must check its result. The CLI exits 1 for partial
growth; it does not erase the partial exports. Fixed-seed replay is tested
within the same executable/environment, not asserted across C libraries,
compiler versions or Ceres versions.

Implicit circles/spheres and cubes now consistently include their boundary,
where the initial root lies. Sampling and root initialization preserve
translated centers; cube membership checks both sides of every axis. These
are domain-contract fixes, not changes to pressure/flow laws or radius tuning.

Local results do not imply green remote CI. The workflow now invokes CTest after
building, but no remote run was submitted as part of this local build repair.
