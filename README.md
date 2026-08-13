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

To generate a single cage with a target number of vertices, run the executable file as

`exeCageGenerator.exe default path-to-input path-to-out-dir target_Nv`

To generate nested cages with target numbers of vertices, run the executable file as

`exeCageGenerator.exe default path-to-input path-to-out-dir target_Nv_0 target_Nv_1 ... target_Nv_n`

The `default` can be substituted by a Phase 2 preset or a configure file. We provide an example at "src/config.json".

Outputs are written under a unique run directory inside the input-name folder:

`path-to-out-dir/input_name/<run_timestamp>__phase2_<mode>[__<mode-specific-details>]/`

The default mode retains its existing `__collapse_hausdorff__flip_<mode>__relocate_<mode>` suffix for output compatibility.

If that directory already exists, the program appends `_001`, `_002`, and so on to avoid overwriting previous results.
