# stochasticConvectionDiffusion-MLE

Source code for the article Maximum likelihood estimation for a stochastic convection-diffusion equation
see [] for details.

The repository encloses code From Matlab, Fortran and Python.

To compile fortran code with the ifx intel compiler and CMake:

- Configure the Build with CMake

```bash

mkdir build && cd build
cmake .. -DCMAKE_Fortran_COMPILER=ifx

```

- Build the Project

```bash

cmake --build .


```

If you build directory is at the same branch of Fotran/src the config with

```bash
cmake ../src -DCMAKE_Fortran_COMPILER=ifx
```
