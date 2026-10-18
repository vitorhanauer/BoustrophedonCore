# Validation record

## Boundary clearance change — 2026-10-18

A new public-API regression first failed because generated waypoints were on
original polygon boundaries. After the change, tests sample the original free
area (including corners and obstacle edges) against the disk footprint at 0%,
25% and 75% overlap. They independently check segment safety and separation
from rectangular obstacles. Narrow passages, domains smaller than the sensor
width, ring winding and fixed boundary terminals are also exercised.

CMake is unavailable in this environment. All five test executables were built
directly with `cc -std=c11 -O2` and the warning flags from CMakeLists.txt; all
passed, including allocator failure injection with and without an obstacle.
The 14-case benchmark matrix and the concave and obstacle fixtures passed.
The largest square case used 3,972 KiB peak RSS; these runs were on the local
host, without a 512 MiB limit or a controlled 1 GHz CPU frequency.

The four functional suites passed AddressSanitizer and UndefinedBehaviorSanitizer
with `ASAN_OPTIONS=detect_leaks=0`. LeakSanitizer could not run under this
environment's ptrace instrumentation; no new leak-check claim is made.

The records below describe the earlier implementation and remain historical.

Validated on 2026-10-15. Implementation source at
`fc25a44bd0bdc20f2548269f73ac5c89d42a32aa`; the following documentation commit
adds these records without changing the implementation.

## TDD and functional checks

Before implementing decomposition, ordering and cell coverage, their suites were
compiled against stubs. CTest reported 0/3 passing: decomposition failed on its
first rectangular fixture, graph on chain connectivity, and sweep on rectangular
coverage. The stitching/API suite was then compiled and failed on its first
`bp_shortest_path` call before the stitching implementation was started.

A further regression reproduced an unsafe diagonal connection when a sweep cut
coincided with a vertical notch. The corrected boundary-arc traversal passed that
test. Inclined polygon fixtures also reproduced rejection of rounded boundary
intersections; coordinate-aware floating-point error bounds corrected it.

Final results:

- Release: **10/10 CTest checks passed** (about 0.17 seconds on the local host).
- ASan, UBSan and LeakSanitizer: **10/10 passed** (6.77 seconds).
- Docker Release build, GCC 12.2: **10/10 passed**.
- Both static and shared libraries built and installed successfully.
- The example links the shared library and includes only the public header.
- All four functional suites, 70 deterministic general polygon fixtures,
  multiple inclined obstacles, a closed start/end route, analytical rectangular
  coverage, percentage spacing, invalid inputs and allocation-failure injection
  are included. Failed allocations leave results empty and leak checks pass.

Reproduce locally:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
ctest --test-dir build --output-on-failure
cmake --install build --prefix /tmp/boustrophedon-install
cmake -S . -B build-sanitize -DCMAKE_BUILD_TYPE=Debug -DBP_SANITIZE=ON
cmake --build build-sanitize -j2
ctest --test-dir build-sanitize --output-on-failure
```

LeakSanitizer was run outside the tool sandbox because it cannot operate under
that sandbox's ptrace instrumentation. No sanitizer was disabled. Sanitizer
memory overhead is excluded from performance measurements. The benchmark's
children use `_exit`; allocation leaks are checked in the ordinary functional
and allocation-failure test executables rather than at benchmark child exit.

## Resource-constrained measurements

The 14 required combinations ran serially, one child process per case, inside
Docker. All returned `success`. Raw results are in [benchmark.csv](benchmark.csv).
Two extra cases, with concavity and an obstacle and 25% overlap, also succeeded;
see [benchmark-extra.csv](benchmark-extra.csv).

| Property | Observed value |
|---|---|
| Image | `sha256:607a533b521998c7785391a650225d8e11e10ac7e742f7004d64b9af038a3ebd` |
| Base image | Debian bookworm-slim |
| Host CPU | AMD Ryzen 5 5500, 12 online logical CPUs |
| Host kernel | Linux 7.1.8-arch1-3, x86_64 |
| Container UID | 10001 |
| `memory.max` | 536870912 bytes (512 MiB) |
| `memory.swap.max` | 0 |
| `cpu.max` | `100000 100000` (one CPU quota) |
| Baseline overlap | 0% |
| Baseline endpoints | `(0,0)` and `(L,L)` |
| Repetitions | One per size/width pair |

The largest required case, 10,000 × 10,000 with width 1, produced **20,006 points**,
used **2,688 KiB peak child RSS**, and took **10.617194 ms** inside the generation
call. This was also the maximum time and RSS in the baseline matrix. These are
single-run observations on this host, not portable timing bounds. Total container
memory includes the runner and runtime, beyond the measured child's RSS; cgroup
limits apply to the whole container.

CPU quota does **not** fix clock frequency. Execution at 1 GHz has not been
established; that requires externally controlled host/VM hardware. No arbitrary
timeout or performance threshold was imposed. The measurements demonstrate
completion within the selected memory/CPU-time quota for these fixtures, not a
resource guarantee for arbitrary vertex counts or near-zero lane spacing.

Reproduce the constrained run:

```sh
./bench/run-container.sh > /tmp/boustrophedon-benchmark.csv
./bench/run-container.sh --size 100 --width 5 --fixture obstacle --overlap 25
```

The script writes the image ID, hardware information, UID and effective cgroup
limits to stderr, and CSV to stdout. The route is heuristic and uses a fixed
sweep direction; these measurements make no claim of globally shortest routes.
