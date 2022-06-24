# SmolyakApprox.jl

This package implements Smolyak's method for approximating multivariate continuous functions.  Two different types of interpolation schemes are allowed: Chebyshev polynomials and piecewise linear.

To install this package you need to type in the REPL

```julia
using Pkg
Pkg.add("SmolyakApprox")
```

Then the package can be used by typing

```julia
using SmolyakApprox
```

Chebyshev polynomials
---------------------

The nodes are computed using Chebyshev-Gauss-Lobatto, with the approximation grid and the multi-index computed by

```julia
grid, multi_ind = smolyak_grid(chebyshev_gauss_lobatto,d,mu,domain)
```

where `d` is the dimension of the function, `mu` is the layer or approximation order, and domain is a 2d-array (2xd) containing the upper and lower bound on each variable.  If domain is not provided, then it is assumed that the variables reside on the [-1,1] interval.  If `mu` is an integer, then an isotropic grid is computed whereas if `mu` is tuple of integers with length `d`, then an anisotropic grid is computed.

With the grid and multi-index in hand, we can compute the weights, or coefficients in the approximation, according to

```julia
weights = smolyak_weights(y,grid,multi_ind,domain)
```

where `y` is a 1d-array containing the evaluations at each grid point of the function being approximated.  Computation of the weights can be made more efficient by computing the inverse interpolation matrix (this generally needs to be done only once, outside any loops)

```julia
inv_interp_mat = smolyak_inverse_interpolation_matrix(grid,multi_ind,domain)
```

with the weights now computed through

```julia
weights = smolyak_weights(y,inv_interp_mat)
```

Lastly, we can evaluate the Smolyak approximation of the function at any point in the domain by

```julia
y_hat = smolyak_evaluate(weights,point,multi_ind,domain)
```

where `point` (a 1d-array) is the point in the domain where the approximation is to be evaluated.

Piecewise linear
----------------

For piecewise linear approximation equidistant nodes are used where the number of nodes is determined according to the Clenshaw-Curtis grid structure: 2^(mu-1)+1

```julia
grid, multi_ind = smolyak_grid(clenshaw_curtis_equidistant,d,mu,domain)
```

Then the weights are computed using

```julia
weights = smolyak_pl_weights(y,grid,multi_ind,domain)
```

and the approximation computed via

```julia
y_hat = smolyak_pl_evaluate(weights,point,grid,multi_ind,domain)
```

Again `mu` can be either an integer or a tuple of integers depending on whether an isotropic or an anisotropic approximation is desired, and the argument `domain` is unnecessary where the grid resides on [-1,1]^d.

Related packages
----------------

- ChebyshevApprox.jl
- HyperbolicCrossApprox.jl
- PiecewiseLinearApprox.jl

References
----------

My primary references when writing this package were:

Judd, K., Maliar, L., Maliar, S., and R. Valero, (2014), "Smolyak Method for Solving Dynamic Economic Models: Lagrange Interpolation, Anisotropic Grid and Adaptive Domain," Journal of Economic Dynamics and Control, 44, pp.92--123.

Klimke, A., and B. Wohlmuth, (2005), "Algorithm 846: spinterp: Piecewise Multilinear Hierarchical Grid Interpolation in MATLAB," ACM Transactions on Mathematical Software, 31, 4, pp.561--579.
