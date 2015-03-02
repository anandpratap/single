# Single - Quasi One Dimensional Finite Volume Solver for Euler Equation

## Features
- Uses backward Newton with exact full jacobian, the jacobian is *efficiently* calculated using code generation by Tapenade
- Adjoint and sensitivity calculations for an arbitary objective function
- OpenMP support, except for the jacobian calculations, which is work in progress

## Building
- Make sure all the path variables in the Makefile are set correctly
- Generate differentiated subroutine using `make adjoint` and then compile using `make`. You might have to create in `bin` directory in the main folder
- Use `make clean` to clean all the object files, and `make clean_adjoint` to clean adjoint related files

## Lessons Learned
- Efficient implementation of a tapenade based solver from scratch, using local stencils to minimize computations
- Nitty gritty details of backward mode of automatic differentation using tapenade
- First substantial FORTRAN program written from scratch
- Makefile stuffs

## Future Enhancements etc.
- OpenMP support for jacobian calculations, difficult because of the way Tapenade handle memory
- Building and Test run on Xeon Co-Phi after step 1
- More modularization and parametrization, remove absolute number from the code
- Adaptive CFL strategies which balances between performance and robustness (CFL of 1000 can be used in some cases while other require as low as 10, specially with shock)
- Applying *best practices*



**Author**: Anand Pratap Singh, anandps@umich.edu

