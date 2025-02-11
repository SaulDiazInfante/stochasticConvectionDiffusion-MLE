# stochasticConvectionDiffusion-MLE

Source code for the article Maximum likelihood estimation for a stochastic convection-diffusion equation
see [] for details.

The repository encloses code From Matlab, Fortran and Python.

To compile fortran code with the ifx intel compiler and CMake:

- Configure the Build with CMake

```

mkdir build && cd build
cmake .. -DCMAKE_Fortran_COMPILER=ifx

```

- Build the Project

```

cmake --build .

```
