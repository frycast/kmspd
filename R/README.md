# kmspd
k-means on SPD matrices

### Instructions

* logchol_v1.R follows the conventions in [lin2019].
* logchol_v2.R recognises that we can leave the diagonals in log space and save a whole bunch of computational hassle. Also, this space contains the 0 matrix.







### References

```
@article{[lin2019],
  title={Riemannian Geometry of Symmetric Positive Definite Matrices via Cholesky Decomposition},
  author={Lin, Zhenhua},
  journal={SIAM Journal on Matrix Analysis and Applications},
  volume={40},
  number={4},
  pages={1353--1370},
  year={2019},
  publisher={SIAM}
}
```