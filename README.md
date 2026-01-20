# RSJM

This repository provides code for simulation and estimation of the **Robust Sparse Jump Model**.

## Main files

- `Utils_sparse_robust_2.R`  
  Contains the main R functions for simulation and estimation.

  Specifically:

- `robust_sparse_jump()`  
  Core function used to estimate the RSJM.

- `robJM.cpp`  
  C++ implementation of computationally intensive routines, integrated into R via **Rcpp**.

  Must be saved in the same folder as Utils_sparse_robust_2.R, and imported via
```r
Rcpp::sourceCpp("robJM.cpp")
