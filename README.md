# SmolyakApprox.jl

This package implements Smolyak's method for approximating multivariate continuous functions.

To install this package (it is currently not official, and so not included in METADATA) you need to run in the REPL

Pkg.clone("https://github.com/RJDennis/SmolyakApprox.jl")

Then the package can be used by typing

using SmolyakApprox

The nodes are computed using either Chebyshev-Gauss-Lobatto or Legendre-Gauss-Lobatto, with the approximation grid and the multi-index computed by the command

grid, multi_ind = smolyak_grid(chebyshev_gauss_lobatto,d,mu,domain)

where d is the dimension of the function, mu is the layer or approximation order, and domain is a 2d-array (2*d) containing the upper and lower bound on each variable.  If domain is not provided, then it is assumed that the variables reside on the [-1,1] interval.  legendre_gauss_lobatto can obviously be used in place of chebyshev_gauss_lobatto.  If mu is an integer, then an isotropic grid is computed whereas if mu is an array of integers with length d, then an anisotropic grid is computed. 

With the grid and multi-index in hand, we van compute the weights, or coefficients in the approximation, according to

weights = smolyak_weights(y,grid,multi_ind,domain)

where y is a 1d-array containing the evaluations at each grid point of the function being approximated.  Computation of the weights can be made more efficients by computing the inverse interpolation matrix

inv_interp_mat = smolyak_inverse_interpolation_matrix(grid,multi_ind,domain)

with the weights now computed through

weights = smolyak_weights(y,inv_interp_mat)

Lastly, we can evaluate the Smolyak approximation of the function at any point in the domain by

y_hat = smolyak_evaluate(weights,point,multi_ind,domain)

where point (a 1d-array) is the point in the domain where the approximation is to be evaluated.

My primary reference when writing this package was:

Judd, K., Maliar, L., Maliar, S., and R. Valero, (2014), "Smolyak Method for Solving Dynamic Economic Models:
Lagrange Interpolation, Anisotropic Grid and Adaptive Domain," Journal of Economic Dynamics and Control, 44, pp.92--123.
