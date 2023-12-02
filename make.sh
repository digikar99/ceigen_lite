
VERSION="v3.4"
SHELL=/data/data/com.termux/files/usr/bin/bash
ARCH=$(uname -m)
unameOut=$(uname)

# Reference: https://stackoverflow.com/questions/394230/how-to-detect-the-os-from-a-bash-script
case "${unameOut}" in
    *Microsoft*)     OS="WSL";; #must be first since Windows subsystem for linux will have Linux in the name too
    *microsoft*)     OS="WSL2";; #WARNING: My v2 uses ubuntu 20.4 at the moment slightly different name may not always work
    Linux*)     OS="linux";;
    Darwin*)    OS="mac";;
    CYGWIN*)    OS="cygwin";;
    MINGW*)     OS="windows";;
    *Msys)     OS="windows";;
    *)          OS="UNKNOWN:${unameOut}"
esac

case "${unameOut}" in
    *Microsoft*)     OSEXT="WSL";; #must be first since Windows subsystem for linux will have Linux in the name too
    *microsoft*)     OSEXT="WSL2";; #WARNING: My v2 uses ubuntu 20.4 at the moment slightly different name may not always work
    Linux*)     OSEXT="linux.so";;
    Darwin*)    OSEXT="mac";;
    CYGWIN*)    OSEXT="cygwin";;
    MINGW*)     OSEXT="windows.dll";;
    *Msys)     OSEXT="windows.dll";;
    *)          OSEXT="UNKNOWN:${unameOut}"
esac

if [ "$ARCH" == "x86_64" ]; then
    if [ "$OS" == "linux" ] || [ "$OS" == "windows"]; then
       # Download the shared library from the release
       URL="https://github.com/digikar99/ceigen_lite/releases/download/$VERSION/libceigen_lite-$ARCH-$OSEXT"
       wget $URL
    else
        export EIGEN_USE_MKL_ALL=1
        g++ --std=c++11 -O3 -mavx2 -shared -o libceigen_lite-$ARCH-$OS.so -fpic ceigen_lite.cpp -I./ -fopenmp
        # gcc -O3 -mavx512f -shared -o libbmas.so -fpic bmas.c
        # gcc -O3 -msse2 -shared -o libbmas.so -fpic bmas.c
        # TODO: Prepare according to different architectures
    fi
else
    ## FIXME for aarch64
    echo "Don't know how to obtain or compile ceigen_lite for architectures other than x86_64. Please issue a PR."
    exit 1
fi
