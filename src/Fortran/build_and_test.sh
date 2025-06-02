#!/bin/bash

set -e

# Set Intel oneAPI environment
source /opt/intel/oneapi/setvars.sh

# Clean build directory
rm -rf build
mkdir build
cd build

# Configure and build
cmake ..
cmake --build .

# Run tests
ctest --output-on-failure
