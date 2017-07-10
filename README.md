# altginv

Matlab ADMM code to compute some norm-minimizing generalized inverses. Written for the following two-part paper series by Ivan Dokmanić and Rémi Gribonval:

- [Beyond Moore-Penrose, Part I: Norm-Minimizing Generalizerd Inverses](https://arxiv.org/abs/1706.08349)
- [Beyond Moore-Penrose, Part II: The Sparse Pseudoinverse](https://arxiv.org/abs/1706.08701)


## Main functions


| Filename          | Description                                           |
| ----------------- | -------------------                                   |
| ginvadmm.m        | (Linearized) ADMM for generalized inverse computation |

## Helper functions

| Filename          | Description                                                   |
| ----------------- | -------------------                                           |
| prox_l1.m         | Proximal operator for the l1 entrywise norm                   |
| prox_l1l1.m       | Proximal operator for the l1->l1 induced norm                 |
| prox_l1l2.m       | Proximal operator for the l1->l2 induced norm                 |
| prox_row21.m      | Proximal operator for the (2,1)-rowwise mixed norm              |
| proj_col21.m      | Projection onto the norm ball of the (2,1)-columnwise mixed norm |

For more details about the relevant definitions see the above two papers. __The proximal mapping for the l1->l1 induced norm uses a projection operator proposed in the following paper:__

> Quattoni, A., Carreras, X., Collins, M., & Darrell, T. (2009). An Efficient Projection for $\ell_{1,\infty}$ Regularization (pp. 857–864). Presented at ICML 2009. http://doi.org/10.1145/1553374.1553484

The file projL1Inf.c is the authors' implementation in C (MEX) and should be MEX-compiled in Matlab on the target platform. 



## Examples and paper figures

| Filename                       | Description                            |
| -----------------              | -------------------                    |
| Example_FrobeniusSpinv.m     | Creates the figures related to the Frobenius norm of the Moore-Penrose pseudoinverse and of the sparse pseudoinverse |
