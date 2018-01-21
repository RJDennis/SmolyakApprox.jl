module SmolyakApprox

include("chebyshev_gauss_lobatto.jl")
include("legendre_gauss_lobatto.jl")
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

export chebyshev_gauss_lobatto,
       legendre_gauss_lobatto,
       smolyak_grid,
       m_i,
       smolyak_weights,
       smolyak_inverse_interpolation_matrix,
       smolyak_evaluate,
       chebyshev_polynomial_derivative,
       smolyak_derivative_finite_difference,
       smolyak_derivative

end
