#!/usr/bin/sh
# Build and install minc-tools in the MouseAtlas tree. Requires libminc
# to be already installed.

MA=/opt/MouseAtlas
#export MA=$HOME/MouseAtlas/Build/
#MA=$HOME/MouseAtlas/Build/debug

export LIBMINC_DIR=$MA/lib

mkdir -p ./build-dir

cmake \
    -B./build-dir \
    -H. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$MA \
    -DCMAKE_COLOR_MAKEFILE=False \
    -DCMAKE_VERBOSE_MAKEFILE=True \
    -DLIBMINC_DIR=$LIBMINC_DIR
