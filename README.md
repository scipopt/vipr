# VIPR: Verifying Integer Programming Results

## About

*VIPR* is a software project to verify, in exact rational arithmetic, the correctness of results computed by mixed-integer linear programming solvers.  It is based on an elementary file format for LP-based branch-and-cut certificates proposed in the article

> Kevin K.H. Cheung, Ambros Gleixner, and Daniel E. Steffy: [Verifying Integer Programming Results](http://dx.doi.org/10.1007/978-3-319-59250-3_13). In: F. Eisenbrand and J. Koenemann, eds., Integer Programming and Combinatorial Optimization: 19th International Conference, IPCO 2017, pp. 148-160, 2017, [`doi:10.1007/978-3-319-59250-3_13`](http://dx.doi.org/10.1007/978-3-319-59250-3_13).

This repository contains a detailed technical specification of the certificate file format in [Version 1.0](cert_spec_v1_0.md) and [Version 1.1](cert_spec_v1_1.md), [software](code/) to check, display, compress, and complete certificate files, and [supplementary information](experiments/) on the computational experiments conducted for the article above.

The specification and handling of incomplete derivations in [Version 1.1](cert_spec_v1_1.md) was added for the purpose of the paper

> Leon Eifler and Ambros Gleixner: [Safe and Verified Gomory Mixed Integer Cuts in a Rational MIP Framework](https://doi.org/10.1137/23M156046X).
 In: SIAM Journal on Optimization 34.1, pp. 742–763.

Please cite both publications if you use VIPR in your work.

## Software

*VIPR* currently provides six C++ scripts, each being called from a terminal together with an appropriate `.vipr` certificate file:

- `vprchck`: A program that verifies mixed-integer linear programming certificate files specified in the `.vipr` file format.
- `vprchck_parallel`: A multi-threaded version of `viprchk`. To ensure highest confidence we continue to support the old single-threaded version.
- `vipr2html`: A program that converts `.vipr` certificate files to a human readable HTML format (not recommended for large files).
- `viprttn`: A program that tightens and improves `.vipr` files, potentially reducing their size and allowing for easier checking.
- `viprcomp`: A program that completes incomplete `.vipr` certificate files in parallel using the exact LP solver `SoPlex`.
- `viprincomp`: A program that makes derivations incomplete. Only useful for testing `viprcomp`.

## File format specification `.vipr`

A conceptual description of the verified integer programming result (`.vipr`) file format is given in the above articles.  A more detailed technical specification is provided [here](cert_spec_v1_1.md).

A small example is given as [paper_eg3.vipr](code/paper_eg3.vipr).  Certificates for large MIP instances from the literature can be found as part of the [supplementary information](experiments/) of the article.

## Installation

The `vipr` scripts are compiled using [CMake](https://cmake.org/).

### CMake Build System

CMake is a build system generator that can create, e.g.,
Makefiles for UNIX and macOS or Visual Studio project files for Windows.

CMake provides an
[extensive documentation](https://cmake.org/cmake/help/latest/manual/cmake.1.html)
explaining available features and use cases as well as an
[FAQ section](https://cmake.org/Wiki/CMake_FAQ). These are the usual steps on a
Linux or macOS system:

    mkdir build
    cd build
    cmake <path/to/vipr>
    make

CMake uses an out-of-source build, i.e., compiled binaries and object files are
separated from the source tree and located in another directory, e.g, `build`.
From within this directory, run `cmake <path/to/vipr>` to configure your build,
followed by `make` to compile the code according to the current configuration.

Afterwards, successive calls to `make` are going to recompile modified source code,
without requiring another call to `cmake`.

Note that in order for `viprcomp` to run, the [SoPlex](https://soplex.zib.de/) and therefore [ZLIB](https://zlib.net/) libraries are required.

If `viprcomp` should not be compiled, it can be turned off in the `cmake <path/to/vipr>` call by using `-DVIPRCOMP=off`.

## How to use VIPR

After installing, run any of the vipr scripts as `./<viprscript> <path/to/.vipr-file>`.

The script `viprcomp` is the only one with the additional option to set verbosity levels as well as the option to disable SoPlex.
The verbosity level of SoPlex can be set to levels 0-5 using the flag `--vebosity=<level>`. Additional debug output can be enabled using `--debugmode=on`.
If it is known that only weak derivations need to be completed, perfomance can be improved by setting `--soplex=off`.

An example call for the completion script: `./viprcomp --verbosity=1 --debugmode=off --soplex=on <path/to/.vipr-file>`.

## Developers and contributors

- [Kevin K.H. Cheung](https://carleton.ca/math/people/kevin-cheung/), School of Mathematics and Statistics, Carleton University
- [Ambros Gleixner](http://www.zib.de/gleixner), Zuse Institute and HTW Berlin
- [Daniel E. Steffy](https://desteffy.github.io), Gurobi Optimization
- [Leon Eifler](https://www.zib.de/members/eifler), Zuse Institute Berlin
- Fabian Frickenstein, TU Berlin

## Software for generating `.vipr` certificates

An [exact rational extension](https://github.com/scipopt/scip/tree/exact-rational) of the solver [SCIP](https://scipopt.org) is able to produce certificates of MIP results.
