# linear-radau
This repository contains fortran90 modules for the solution of stiff, linear, first order initial value problems. The routines are based on a linearised version of the RadauIIA described in [Hairer & Wanner](http://link.springer.com/book/10.1007%2F978-3-642-05221-7). They numerically solve the initial value problem y' = Ay with given initial conditions.

## Dependencies
- [MUMPS] (http://mumps.enseeiht.fr/)

## Compilation
An example makefile is provided. You'll first need to compile MUMPS and link the radau code to your MUMPS installation afterwards

## Usage
All versions of the radau integrator use the assembled data format to provide the matrix to MUMPS. This format uses three arrays: one for the nonzero values of the matrix, one for the row indices and one for the column indices. Some basic matrix I/O routines are also provided.
