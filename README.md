# sinc-colloc-vol-fred
Numerical solvers for Volterra-Fredholm integral equations of the second kind by Sinc-collocation methods

## Overview
These programs solve three examples of Volterra-Fredholm integral equations
of the second kind, conducted in [2].

Those problems are solved by means of the following 4 methods:
* Original SE-Sinc-collocation method [3]
* Original DE-Sinc-collocation method [1]
* New SE-Sinc-collocation method [2]
* New DE-Sinc-collocation method [2]

The name of the program denotes the method and example number. For
example, SE_orig_ex1.c denotes the Original SE-Sinc-collocation method
for Example 1, and DE_new_ex2.c denotes the New SE-Sinc-collocation method
for Example 2.

LAPACK in Apple's Accelerate framework is used for computation of the
system of linear equations and its condition number. If you want to use
another LAPACK library, modify make files according to your installation.

Each program solves those problems increasing N as N = 5, 10, 15, 20, ...,
and outputs n=2N+1, maximum error over the target interval, and the condition
number (with the infinity norm) of the resulting system.

## Results
Outputs by those programs are stored in data/ directory, with .dat extension.
Gnuplot programs for creating graphs are also stored in the directory.

computation environment:

OS: macOS Monterey  
CPU: 2.4 GHz (quad core) Intel Core i5  
Memory: 16 GB 2133 MHz LPDDR3  
Compiler: Apple clang version 14.0.0  
Library: LAPACK (Apple's Accelerate framework)

## References
1. E. D. John and N. Ogbonna: A double exponential Sinc collocation method
 for Volterra-Fredholm integral equations of the second kind, J. Nigerian
 Math. Soc., Vol. 35 (2016), pp. 408--423.
2. T. Okayama: Improvement of Sinc-collocation methods for Volterra-Fredholm
 integral equations of the second kind and their theoretical analysis,
 arXiv:2503.11569 [math.NA], Mar 2025.
3. A. S. Shamloo, S. Shahkar, and A. Madadi: Numerical solution of the
 Fredholme-Volterra integral equation by the Sinc function, Am. J. Comput.
 Math., Vol. 2 (2012), pp. 136--142.
