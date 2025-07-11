#!/bin/bash

# Set script to exit on error
set -e

NCORES=4
unamestr=`uname`
if [[ "$unamestr" == "Linux" ]]; then
        NCORES=`grep -c ^processor /proc/cpuinfo`
fi

if [[ "$unamestr" == "Darwin" ]]; then
        NCORES=`sysctl -n hw.ncpu`
fi

# Define build and deploy directory
BUILD_DIR="build"
DEPLOY_DIR="deploy"

# Remove old build directory if it exists
if [ -d "$BUILD_DIR" ]; then
    rm -rf "$BUILD_DIR"
fi

# Remove old deploy directory if it exists
if [ -d "$DEPLOY_DIR" ]; then
    rm -rf "$DEPLOY_DIR"
fi

# Create a new build directory
mkdir "$BUILD_DIR"
cd "$BUILD_DIR"

# Run CMake and Build
cmake ..

cmake --build . --parallel $NCORES

cd ..

# Create a new deploy directory
mkdir "$DEPLOY_DIR"

cp ./build/multi_machine_scheduler deploy/multi_machine_scheduler
