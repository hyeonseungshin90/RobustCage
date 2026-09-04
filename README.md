# README

This is the reference implementation of "Robust Coarse Cage Construction with Small Approximation Errors".

## Dependencies

* Boost 1.79 (The only library you need to install on your computer)
* Other libraries are contained src/ThirdPartyLib

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

Outputs are written under a unique run directory inside the input-name folder:

`path-to-out-dir/input_name/<run_timestamp>[__phase1_topological_offset]__phase2_<mode>[__<mode-specific-details>]/`

The default mode retains its existing `__collapse_hausdorff__flip_<mode>__relocate_<mode>` suffix for output compatibility.

If that directory already exists, the program appends `_001`, `_002`, and so on to avoid overwriting previous results.
