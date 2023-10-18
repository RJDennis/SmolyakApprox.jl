# SmolyakApprox.jl

Introduction
============

This package implements Smolyak's method for approximating multivariate continuous functions.  Two different types of interpolation schemes are allowed: Chebyshev polynomials and piecewise linear. The package also implements Clenshaw-Curtis integration and Gauss-Chebyshev quadrature.

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

The nodes are computed using Chebyshev-Gauss-Lobatto (Chebyshev extrema), with the approximation grid and the multi-index computed by

```julia
grid, multi_ind = smolyak_grid(chebyshev_gauss_lobatto,d,mu,domain)
```

where `d` is the dimension of the function, `mu` is the layer or approximation order, and domain is a 2d-array (2xd) containing the upper and lower bound on each variable.  If domain is not provided, then it is assumed that the variables reside on the [-1,1] interval.  If `mu` is an integer, then an isotropic grid is computed whereas if `mu` is a 1d array of integers with length `d`, then an anisotropic grid is computed.  Because the chebyshev_gauss_lobatto points are the same as the Chebyshev extrema points you can use `chebyshev_extrema` in place of `chebyshev_gauss_lobatto`. 

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

You can evaluate the Smolyak approximation of the function at any point in the domain by

```julia
y_hat = smolyak_evaluate(weights,point,multi_ind,domain)
```

where `point` (a 1d-array) is the point in the domain where the approximation is to be evaluated.

Lastly, you can compute derivatives, gradients, and hessians according to:

```julia
d = smolyak_derivative(weights,point,multi_ind,domain,pos) # Takes the derivative with respect to variable 'pos'
g = smolyak_gradient(weights,point,multi_ind,domain)
h = smolyak_hessian(weights,point,multi_ind,domain)
```

Integration
-----------

To numerically integrate a function, you first create an approximation plan and then call the integration function:

```julia
plan = smolyak_plan(chebyshev_gauss_lobatto,d,mu,domain)

integral = smolyak_integrate(f,plan,:clenshaw_curtis)      # uses Clenshaw-Curtis
integral = smolyak_integrate(f,plan,:gauss_chebyshev_quad) # uses Gauss-Chebyshev quadrature
```
where `f` is the function to be integrated and `plan` is the approximation plan, discussed below.  Both methods integrate the function over the full approximation domain.

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

Again `mu` can be either an integer or a 1d array of integers depending on whether an isotropic or an anisotropic approximation is desired, and the argument `domain` is unnecessary where the grid resides on [-1,1]^d.

Multi-threading
---------------

There are multi-threaded functions to compute the polynomial weights and the interpolation matrix.  These multi-threaded functions are accessed by adding `_threaded` to the end of the funtion, as per

```julia
weights = smolyak_weights_threaded(y,inv_interp_mat)
```

Useful structures
-----------------

The key structure to be aware of is the SApproxPlan, which contains the key information needed to approximate a function.

```julia
d = 3
mu = 3
domain = [2.0 2.0 2.0; -2.0 -2.0 -2.0]
grid, mi = smolyak_grid(chebyshev_extrema,d,mu,domain)
plan = SApproxPlan(:chebyshev_extrema,grid,mi,domain)
```
or
```julia
plan = smolyak_plan(chebyshev_extrema,d,mu,domain)
```

Once the approximation plan has been constructed it can be used to create functions to interpolate and to compute gradients and hessians.

```julia
f = smolyak_interp(y,plan)
g = smolyak_gradient(y,plan)
h = smolyak_hessian(y,plan)

point = [1.0, 1.0, 1.0]

f(point)
g(point)
h(point)
```

There are threaded versions of `smolyak_interp`, `smolyak_gradient`, and `smolyak_hessian`; just add `_threaded` to the end of the function name.

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
