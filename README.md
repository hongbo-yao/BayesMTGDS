# BayesMTGDS

BayesMTGDS is a standalone C++ code for trans-dimensional Bayesian joint inversion of magnetotelluric (MT) and geomagnetic depth sounding (GDS) responses. Supported data and transfer functions include:

- MT apparent resistivity and phase

- GDS global Qn and Cn responses (a global network of geomagnetic observatory data, satellite data; magnetospheric ring current and ionospheric Sq current (not recommended))

- GDS local Cn responses (single site ground-based data; magnetospheric ring current and ionospheric Sq current (not recommended))

- GDS global-to-local transfer functions (single site ground-based data; magnetospheric ring current and ionospheric Sq current)

### Reference
The background theory is described in:

- Hongbo Yao, Zhengyong Ren*, Jingtian Tang, Rongwen Guo, Jiayong Yan. Trans-dimensional Bayesian joint inversion of magnetotelluric and geomagnetic depth sounding responses to constrain mantle electrical discontinuities. Geophysical Journal International, 2023, 233(3), 1821-1846. https://doi.org/10.1093/gji/ggad029

If this helps you, a citation to our paper is appreciated.

### Authors

- Hongbo Yao, Central South University, Email: yaohongbo@csu.edu.cn, https://www.researchgate.net/profile/Hongbo-Yao-2

- Zhengyong Ren, Central South University, Email: renzhengyong@csu.edu.cn

### Usage
Install MPI and type 'make' to compile the source code. Then please see /examples for details.
