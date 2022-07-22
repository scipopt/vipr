# ReadMe

# VIPR: Verifying Integer Programming Results
---
## About
---
*VIPR* is new a software project to verify, in exact rational arithmetic, the correctness of results computed by mixed-integer linear programming solvers.  It is based on an elementary file format for LP-based branch-and-cut certificates proposed in the article

> Kevin K.H. Cheung, Ambros Gleixner, and Daniel E. Steffy: [Verifying Integer Programming Results](http://dx.doi.org/10.1007/978-3-319-59250-3_13). In: F. Eisenbrand and J. Koenemann, eds., Integer Programming and Combinatorial Optimization: 19th International Conference, IPCO 2017, pp. 148-160, 2017, [`doi:10.1007/978-3-319-59250-3_13`](http://dx.doi.org/10.1007/978-3-319-59250-3_13).

This repository contains a detailed technical [specification of the certificate file format](cert_spec_v1_0.html), [software](code/) to check, display, and compress certificate files, and [supplementary information](experiments/) on the computational experiments conducted for the article above.

## Software
---
*VIPR* currently consists of five C++ scripts, each being called from a terminal together with an appropriate `.vipr` [certificate file](cert_spec_v1_0.html):

- `vprchck`: A program that verifies mixed-integer linear programming certificate files specified in the `.vipr` file format.
- `vipr2html`: A program that converts `.vipr` certificate files to a human readable HTML format (not recommended for large files).
- `viprttn`: A program that tightens and improves `.vipr` files, potentially reducing their size and allowing for easier checking.
- `viprwidth`: A program that returns the cutwidth on derived constraints.
- `viprcomp`: A program that completes incomplete `.vipr` certificate files using `SoPlex`.

## File format specification `.vipr`
---
A conceptual description of the verified integer programming result (`.vipr`) file format is given in the article [above](http://nbn-resolving.de/urn:nbn:de:0297-zib-61044).  A more detailed technical specification is provided [here](cert_spec_v1_0.html).

A small example is given as [paper_eg3.vipr](code/paper_eg3.vipr).  Certificates for large MIP instances from the literature can be found as part of the [supplementary information](experiments/) of the article.
## Installation
---

The `vipr`-scripts are compiled using [CMake](https://cmake.org/).

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

If it is not desired to compile `viprcomp`, it can be turned off in the `cmake <path/to/vipr>` call by using `-DVIPRCOMP=off`.

## How to use
---
After installing, run any of the vipr-scripts by running `./<viprscript> <path/to/.vipr-file>`.

### Verbosity settings for `viprcomp`

`viprcomp` allows for additional verbosity settings for SoPlex as well as the script itself at runtime.
By default, both verbosity settings are turned off.
To change this use the command `./viprcomp <SoPlex verbosity level> <Debugmode ON/OFF> <path/to/.viprfile>`.
Both settings are optional and may be left out.

Soplex options are

    0   Error (default)
    1   Warning
    2   Debug
    3   Normal
    4   High
    5   Full

For further information, consult the [SoPlex Documentation](https://soplex.zib.de/doc/html/).

Debugmode is either set to `ON` or `OFF` (default).


## Developers
---
- [Kevin K.H. Cheung](https://carleton.ca/math/people/kevin-cheung/), School of Mathematics and Statistics, Carleton University
- [Ambros Gleixner](http://www.zib.de/gleixner), Department Optimization, Zuse Institute Berlin
- [Daniel E. Steffy](https://files.oakland.edu/users/steffy/web/), Department of Mathematics and Statistics, Oakland University
- [Fabian Frickenstein](https://www.zib.de/members/frickenstein), Mathematical Optimization Methods, Zuse Institute Berlin

## Software for generating `.vipr` certificates
---
We have created an extension of the [exact rational version](http://scip.zib.de/#exact) of the [SCIP optimization software](http://scip.zib.de) that is able to produce certificates of MIP results.