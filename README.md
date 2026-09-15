# README

This is the reference implementation of "Robust Coarse Cage Construction with Small Approximation Errors".

## Dependencies

* Boost 1.79 (The only library you need to install on your computer)
* Other libraries are contained src/ThirdPartyLib
* The bundled Geometry Central v1.1.0 supplies intrinsic flip geodesics for
  initial boundary-rail construction.

## How to build

Configure the CMake file and build it.

We succeed to build our code on Windows10 by MSVC-16-x64.

We also provide a binary executable file in supplementary.

## How to use

Input meshes may be supplied as OBJ, OFF, PLY, OM, or STL files. Both ASCII
and binary STL are supported, including binary STL files whose 80-byte header
starts with `solid`.

To generate a single cage with a target number of vertices, run the executable file as

`exeCageGenerator.exe default path-to-input path-to-out-dir target_Nv`

For example, to generate a cage directly from an STL file:

`exeCageGenerator.exe default path-to-input.stl path-to-out-dir target_Nv`

To generate nested cages with target numbers of vertices, run the executable file as

`exeCageGenerator.exe default path-to-input path-to-out-dir target_Nv_0 target_Nv_1 ... target_Nv_n`

The `default` can be substituted by a Phase 2 preset or a configure file. We provide an example at "src/config.json".

Phase 1 has two modes:

* `default` keeps the original two rounds of global 1-to-12 tetrahedron subdivision.
* `phase1_topological_offset` uses the simplicial embedding and offset insertion of Zint et al., then extracts the initial-cage boundary.

For example:

`exeCageGenerator.exe phase1_topological_offset path-to-input path-to-out-dir target_Nv`

Phase presets compose with `+`, for example:

`exeCageGenerator.exe phase1_topological_offset+phase2_linear_solve path-to-input path-to-out-dir target_Nv`

Linear-solve Phase 2 repeats `collapse -> flip -> rail update` for up to 30 cycles, then runs
one final relocation stage. Each collapse stage rebuilds its candidates after
the previous cycle's flips and rail relabeling, so newly feasible collapses can be accepted.
The global target edge length remains fixed across cycles. Once the vertex
target is reached, collapse is skipped while flips and rail updates can continue.
Rail updates run only when boundary rails are enabled. A cycle with no accepted
collapses, flips, or rail updates ends this loop early. Relocation starts after
the loop ends; no collapse, flip, or rail update follows it.
Collapse placement and its line-search energy use plane approximation and the triangle-quality
surrogate. Let `L` be the local mean edge length for a collapse candidate.
The plane term sums the area-weighted squared plane distances and divides by `L^4`:
`E_plane(x) = planeWeight * sum_f A_f * d_f(x)^2 / L^4`.
The triangle surrogate averages the squared Euclidean distances to the
equilateral reference apices and divides by `L^2`:
`E_triangle(x) = triangleQualityWeight / (N * L^2) * sum_f ||x - reference_apex_f||^2`.
Each triangle has equal weight in the surrogate, with equal penalties in every
direction. Both terms are dimensionless, so uniform scaling of the geometry
preserves their values, apart from numerical safeguards. Their relative strengths
are set by `planeWeight` and `triangleQualityWeight`. Older JSON files using
`qemWeight` are accepted as a fallback; `planeWeight` takes precedence when
both are present. Serialized settings use `planeWeight`. Placement, line search and queue
scoring use the same normalized terms. Uniformity contributes only to collapse
queue priority: the global mode adds
`uniformityWeight * (edge_length / target_length)^2` to the score.
Flips prioritize improvement in the minimum quality of their two triangles
and reject intersections, as described below.
Relocation combines a tangential smoothing target `t`, computed
from neighboring vertices' mixed Voronoi areas (Botsch-Kobbelt area-equalizing
smoothing), with the source surface's closest point `s` to the current vertex.
Both targets remain fixed during that vertex's line search. The target is their
weighted average, minimizing
`E(x) = tangentialWeight * ||x-t||^2 + surfaceWeight * ||x-s||^2`.
Backtracking from the current position toward this target uses
`alpha = 1, 1/2, 1/4, ...`. An accepted move must strictly decrease this energy,
pass the degeneracy, orientation and intersection checks, and preserve
`new_min_quality >= min(old_min_quality, minTriangleQuality)` within roundoff
tolerance. Quality is normalized to `[0, 1]`; its default floor is `0.2`.
A fan above the floor may decrease to the floor; a fan below it cannot worsen.
With `+boundary_rail`, a mixed rail/nonrail edge may collapse only from the
nonrail vertex into the existing rail vertex at its unchanged position. This
preserves the rail label and loop connectivity and must pass the existing
collapse validity checks. Actual rail-edge collapses still project their new
position onto the source-boundary support half-strips. Rail edges cannot flip.
Rail vertices and actual mesh boundary vertices remain fixed during relocation.
Flips and relocation preserve vertex count.

The linear-solve flip stage uses triangle quality for priority and acceptance,
including when boundary rails are enabled. Every accepted flip requires a
finite increase greater than `1e-12` in the minimum quality of its two triangles.
Rail chord counts do not affect priority or block a flip. All candidates
retain the existing topology, valence, degeneracy and source/cage intersection
checks; there is no additional flip quality floor. Changed neighborhoods are
requeued and candidates are rechecked before acceptance. Flips preserve rail
labels, rail vertex positions, rail edges and rail perimeter, but may change
the cage surface. Actual rail-edge collapse, its projection and line search
remain unchanged.

The C++ entry point is `do_quality_flip()`, with no chord-policy argument or
optional chord-aware mode. The triangle-based rail update runs after flipping
when boundary rails are enabled.

Initial boundary rails are constructed before Phase 2 with
`Dijkstra -> intrinsic flip geodesics -> actual mesh splitting`.
After inserting the source-boundary anchors, the builder finds closed, simple
edge paths with Dijkstra and its existing path validation. Geometry Central
then shortens the rail network intrinsically, keeping every anchor fixed.
The shortened paths are traced onto the cage surface; crossing points become
actual cage vertices, and crossed edges and faces are split so each rail is a
chain of real mesh edges usable by the existing collapse implementation.
Embedding keeps existing vertex positions and the cage surface fixed, but
can add vertices and triangles. It does not replace a path across folded faces
with a straight 3D chord.

This construction runs for both `boundary_rail` and `boundary_rail_vertex`,
including the independent builds in `boundary_rail_compare`, without a new
command-line option. If intrinsic shortening or embedding fails, the builder
keeps the cage and rails from the completed Dijkstra construction and records
a warning and failure statistics; that result has not received the geodesic
update. Flip geodesics seek locally shortest paths, not a global optimum, and
fixed anchors can leave a path unchanged. This step is performed only during
initial rail construction. Later collapses do not preserve a geodesic
guarantee or trigger another intrinsic shortening pass.

With boundary rails enabled, linear-solve Phase 2 calls the triangle-based
`update_boundary_rails()` after each quality-flip pass. Each call repeatedly
replaces a labeled path `A-B-C` with `A-C` across a cage triangle and releases
`B`, until no eligible shortcut remains. It preserves a simple loop of at least
three vertices, changes labels only, and has no source-distance or opening-size
error bound. Relabeling alone counts as progress, so the next cycle can use the
new rail constraints even when the preceding collapse and flip counts were zero.

The final relocation stage performs up to 20 full-mesh sweeps, recomputing the
Voronoi areas, normals and source closest points from current positions. A sweep
accepting no moves ends relocation early. JSON settings control these three
separate limits:

* `paramCageSimplifier.phase2QualityPolishIterations`: collapse/flip/rail-update cycles,
  default `30`. The existing key name is retained; `0` runs collapse once and
  disables flips, rail updates, and final relocation.
* `paramCageSimplifier.paramRelocate.qualitySweeps`: final relocation sweep limit,
  default `20`; `0` disables relocation.
* `paramCageSimplifier.paramRelocate.lineSearchMaxIter`: backtracking attempts
  per vertex per sweep, default `12`; `0` also disables relocation.

The relocation energy and quality gate have three further JSON settings under
`paramCageSimplifier.paramRelocate`:

* `tangentialWeight`, default `1.0`, and `surfaceWeight`, default `1.0`: finite,
  nonnegative weights. Zero surface weight gives the tangential target only;
  zero tangential weight gives the source closest-point target only. Both zero
  disables relocation.
* `minTriangleQuality`, default `0.2`: finite quality floor in `[0, 1]`.

The sweep cap allows repeated smoothing and surface attraction while early
stopping avoids spending all 20 sweeps on a stationary mesh. The local energy
does not minimize or guarantee improvement of the full Hausdorff distance.
Geometry checks evaluate each candidate configuration; they do not perform
continuous collision detection or enforce a positive surface clearance.

Outputs are written under a unique run directory inside the input-name folder:

`path-to-out-dir/input_name/<run_timestamp>[__phase1_topological_offset]__phase2_<mode>[__<mode-specific-details>]/`

The default mode retains its existing `__collapse_hausdorff__flip_<mode>__relocate_<mode>` suffix for output compatibility.

If that directory already exists, the program appends `_001`, `_002`, and so on to avoid overwriting previous results.

With `+boundary_rail`, the rails are exported next to each cage OBJ in two
states:

* `input_name_cage_<label>_initial_rails.obj` / `.txt`: the constructed rails,
  including successful geodesic embedding, written before Phase 2 starts.
* `input_name_cage_<label>_rails.obj` / `.txt`: the rails carried by the final
  cage.

Each OBJ holds the rail vertices as OBJ points and the rail edges as OBJ line
elements, one object/group per rail id. Each TXT lists the same rails with the
1-based vertex indices of the mesh named in its header: the final rails index
`input_name_cage_<label>.obj`, while the initial rails index the cage after
rail construction and before Phase 2, so the two index spaces differ.
Rail ids are shared between both states.

The quality and boundary-rail geodesic regression tests can be built and run with:

```powershell
cmake -S src -B src/build -DCAGE_BUILD_TESTS=ON "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
cmake --build src/build --config Release --target exeCageGenerator cageQualityPolishTests cageBoundaryRailGeodesicTests
ctest --test-dir src/build -C Release --output-on-failure
```
