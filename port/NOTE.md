# Notes on Citlali Port


## Benchmarks

### Model fit

The ceres-solver is used for model fit.

7000~ fitting of I/Q data for resonance frequencies take less than 1 second.


### PCA clean

The Eigen/SelfAdjointEigenSolver takes 1s for a 1000x1000 matrix without multi-processing.
The algorithm scales as N^3.

Need to use a better parallelized library/algorithm to speed up.

### 1D FFT

There is FFT involved in the sensitivity calculation.
It is very fast. TODO: add an accurate number.

