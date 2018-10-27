module SmolyakApprox
using LinearAlgebra

include("chebyshev_gauss_lobatto.jl")
include("legendre_gauss_lobatto.jl")
include("clenshaw_curtis_equidistant.jl")
include("chebyshev_polynomial.jl")
include("smolyak_weights.jl")
include("smolyak_evaluate.jl")
include("smolyak_grid.jl")
include("generate_multi_index.jl")
include("combine_nodes.jl")
include("m_i.jl")
include("chebyshev_polynomial_derivative.jl")
include("smolyak_derivative_finite_difference.jl")
include("smolyak_derivative.jl")
include("scale_nodes.jl")
include("smolyak_full_grid.jl")
include("smolyak_piecewise_linear.jl.jl")

export chebyshev_gauss_lobatto,
       legendre_gauss_lobatto,
       clenshaw_curtis_equidistant,
       smolyak_grid,
       m_i,
       smolyak_weights,
       smolyak_inverse_interpolation_matrix,
       smolyak_evaluate,
       chebyshev_polynomial_derivative,
       smolyak_derivative_finite_difference,
       smolyak_derivative,
       smolyak_grid_full,
       smolyak_weights_full,
       smolyak_evaluate_full,
       smolyak_derivative_full,
       smolyak_pl_weights,
       smolyak_pl_evaluate

end
