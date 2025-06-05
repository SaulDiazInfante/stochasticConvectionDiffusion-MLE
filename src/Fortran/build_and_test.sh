#!/bin/bash

set -e

# Clean build directory
rm -rf build
mkdir build
cd build

# Configure and build
cmake ..
cmake --build .

# Run tests
ctest --output-on-failure
