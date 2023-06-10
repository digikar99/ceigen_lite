
SHELL=/data/data/com.termux/files/usr/bin/bash
ARCH=$(uname -m)

if [ "$ARCH" == "x86_64" ]; then
    export EIGEN_USE_MKL_ALL=1
    g++ --std=c++11 -O3 -mavx2 -shared -o libceigen_lite.so -fpic ceigen_lite.cpp -I./ -fopenmp
    # gcc -O3 -mavx512f -shared -o libbmas.so -fpic bmas.c
    # gcc -O3 -msse2 -shared -o libbmas.so -fpic bmas.c
    # TODO: Prepare according to different architectures
else
    ## FIXME for aarch64
    gcc -O3 -shared -o libbmas.so -fpic bmas.c
fi
