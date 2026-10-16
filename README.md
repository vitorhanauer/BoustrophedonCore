# Boustrophedon core

C11 coverage-path library for a simple polygon with disjoint polygonal obstacles.
The public `generateWaypoints` function returns a sequence of Cartesian points
starting and ending at the caller's coordinates. Both static and shared libraries
are built; the header can also be included from C++.

## Build and test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
ctest --test-dir build --output-on-failure
./build/bp_example
cmake --install build --prefix /tmp/boustrophedon-install
```

See `include/boustrophedon.h` and `examples/basic.c` for the C interface. Link with
`-lboustrophedon -lm`. Other languages can bind the same C ABI: point coordinates
are `double`, counts use `size_t`, polygons borrow their input arrays, and returned
paths own their point arrays. No wrapper or coordinate-system conversion is
included. All coordinates and widths must use the same Cartesian length unit.

Initialize the output to `{0}`. On success, call `freeWaypoints(&output)` before
reusing it; freeing an empty result is safe and resets its pointer and count.
Failures return a `BpStatus` and leave output empty. Do not pass a live result as
output without freeing it first. Inputs remain unchanged, and calls do not use
shared mutable state.

## Geometry and coverage model

- Outer and obstacle rings must be simple and have positive area. Either winding
  is accepted, as is a repeated closing vertex. Consecutive duplicate vertices,
  self-intersections, overlapping/nested obstacles and touching contours are
  rejected. Collinear vertices without backtracking are accepted.
- Obstacles must be strictly inside the outer ring. Start and end must be in the
  free domain. Outer and obstacle boundaries are traversable; obstacle interiors
  and the exterior of the area are forbidden. The agent has no physical radius
  or minimum turning radius.
- The coverage footprint is a disk of radius `width/2`, with no sensor occlusion.
  Width is the full diameter. Overlap is a percentage in `[0,100)`, and lane
  spacing is `width * (1 - overlap_percent/100)`. Exactly 100% has zero spacing
  and is rejected. Nonfinite values and unrepresentable spacing are rejected.
- Decomposition and coverage use the vertical axis for the sweep segments and
  advance along x. There is no search over sweep angles in this release.

The decomposition intersects the free domain between vertex x-events and joins
fragments until a split or merge in connectivity. Cells keep their actual
polygonal boundaries, including concavities. Adjacent cells share a boundary
interval, rather than only a single point.

Within each x-monotone cell, the path traverses complete vertical cuts with gaps
no greater than the requested spacing and links successive cuts along the cell
boundary. After the last cut, it connects to the next cell or the final point
without adding a full contour lap. Necessary boundary links remain to avoid
shortcuts outside concave cells.

Full coverage by a disk of radius `width/2` is verified for rectangular cells.
For arbitrary inclined or sharp boundaries, lane spacing alone does not guarantee
that every boundary point is within that radius of the route. No complete
boundary-coverage guarantee is made for those geometries.

The planner starts with nearest-neighbor cell ordering, chooses among four
coverage variants by dynamic programming, and applies at most ten 2-opt passes
that reduce total length. Costs include coverage, safe connectors, and the fixed
start/end. It is deterministic and seeks a short route; it does not guarantee a
global minimum. Cells are covered once, but connectors can revisit them.
Connectors use Euclidean shortest paths in a polygonal visibility graph.

## Numerical and resource limits

Geometry uses `double` coordinates with `long double` intermediates for predicates
and intersections. This is floating-point geometry, not an exact-arithmetic
kernel. Degenerate or numerically unresolvable inputs may return a geometry or
representation error. Recenter and rescale data when coordinate magnitudes are
large compared with feature sizes. There is no grid allocation proportional to
area or a fixed maximum polygon dimension.

Memory grows with vertices, decomposition fragments and generated waypoints.
Visibility queries do not create a matrix involving all lane endpoints. Very
large vertex counts and extremely small spacing can still require substantial
CPU or memory; allocation and representation failures are reported. No runtime
bound is claimed for arbitrary numbers of vertices or obstacles.

The reference is [Coverage Path Planning: The Boustrophedon Decomposition,
Choset and Pignon](https://publications.ri.cmu.edu/coverage-path-planning-the-boustrophedon-decomposition).
See `docs/validation.md` for executed validation and resource measurements.

## Benchmarks and Docker

```sh
./build/bp_runner --all
./build/bp_runner --size 100 --width 5 --fixture obstacle --overlap 25 --repeat 3
./bench/run-container.sh > /tmp/boustrophedon-benchmark.csv
```

The default runner executes all 14 size/width combinations from `AGENTS.md`.
`--fixture concave` and `--fixture obstacle` select additional geometries;
`--help` lists all options. Each case and repetition runs in a fresh child process.
CSV includes elapsed generation time, waypoint count, total Euclidean path
length and the child's peak RSS in KiB (Linux `getrusage`). The timer excludes
fixture setup, process startup and CSV formatting. There is no imposed timeout.

`bench/run-container.sh` builds the image, writes hardware/image metadata to
stderr and CSV to stdout. It sets 512 MiB RAM, disables additional swap and limits
CPU time to one logical CPU. An accessible Docker daemon is required. You can
also run the image directly:

```sh
docker build -t boustrophedon-bench .
docker run --rm --memory=512m --memory-swap=512m --cpus=1 boustrophedon-bench --all
```

The image runs as UID 10001. Docker CPU quotas limit processing time rather than
clock frequency: these measurements do **not** establish performance at 1 GHz.
See [Docker resource constraints](https://docs.docker.com/engine/containers/resource_constraints/).
Use a host or VM with externally controlled frequency for that additional test.

For address, undefined-behavior and leak checks:

```sh
cmake -S . -B build-sanitize -DCMAKE_BUILD_TYPE=Debug -DBP_SANITIZE=ON
cmake --build build-sanitize -j2
ctest --test-dir build-sanitize --output-on-failure
```

LeakSanitizer needs an environment without a ptrace-based sandbox. The allocator
fault-injection suite uses GNU-style linker wrapping and is enabled on Linux with
GCC/Clang. The library itself uses only C11 and the math library; the benchmark
runner uses POSIX process and timing APIs.
