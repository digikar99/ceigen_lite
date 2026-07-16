See [./ceigen_lite.h](./ceigen_lite.h) for a list of available functions.

The argument pattern is similar to BLAS/LAPACK: the first few arguments correspond to the size of the vectors/matrices, the vector/matrices themselves are passed in as C pointers. Some functions also take a layout argument, which should be the character 'R' or 'r' if the pointer denotes a matrix with a row major layout, and should be the character 'C' or 'c' if the pointer denotes a column major layout.

As of this writing, [Eigen](https://gitlab.com/libeigen/eigen/-/tree/master) version 3.4 is included, [licensed under MPL2](https://gitlab.com/libeigen/eigen/-/blob/master/COPYING.README).

This repository also packs the [EigenRand](https://github.com/bab2min/EigenRand), version 0.5.1, licensed under MIT, and it was in fact the primary motivation for this project. EigenRand offers several sampling functions optimized using SIMD. These include sampling from a gaussian/normal distribution, beta distribution, chi-squared, and several others.

You can grab a binary from the [Releases](https://github.com/digikar99/ceigen_lite/releases).

Run `make` (or `gmake`) in the current directory to build the `libceigen_lite-*` shared library. The  full name is platform specific.
