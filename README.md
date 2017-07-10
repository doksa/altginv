# altginv

Matlab ADMM implementation of some norm-minimizing generalized inverses. It was written for the following two-part paper series by Ivan Dokmanić and Rémi Gribonval:

- Beyond Moore-Penrose: 

This repository contains a collection of Matlab routines to accompany the paper
_Euclidean Distance Matrices: A Short Walk Through Theory, Algorithms and Applications_ by Ivan
Dokmanic, Reza Parhizkar, Juri Ranieri and Martin Vetterli.

## Main functions


| Filename          | Description                                           |
| ----------------- | -------------------                                   |
| ginvadmm.m        | (Linearized) ADMM for generalized inverse computation |

## Helper functions

The proximal operator for the l1->l1 induced norm uses a projection operator proposed in the following paper:

Quattoni, A., Carreras, X., Collins, M., & Darrell, T. (2009). An Efficient Projection for $\ell_{1,\infty}$ Regularization (pp. 857–864). Presented at ICML 2009. http://doi.org/10.1145/1553374.1553484

The file projL1Inf.c is the authors' implementation in C (MEX) and should be MEX-compiled in Matlab on the target platform. 


| Filename          | Description                                                   |
| ----------------- | -------------------                                           |
| prox_l1.m         | Proximal operator for the l1 entrywise norm                   |
| prox_l1l1.m       | Proximal operator for the l1->l1 induced norm                 |
| prox_l1l2.m       | Proximal operator for the l1->l2 induced norm                 |
| prox_row21.m      | Proximal operator for the (2,1)-rowwise mixed norm              |
| proj_col21.m      | Projection onto the norm ball of the (2,1)-columnwise mixed norm |

For more details about the relevant definitions see the papers.

## Examples and paper figures

| Filename                       | Description                            |
| -----------------              | -------------------                    |
| Example_FrobeniusSpinv.m     | Creates the figures related to the Frobenius norm of the Moore-Penrose pseudoinverse and of the sparse pseudoinverse |
