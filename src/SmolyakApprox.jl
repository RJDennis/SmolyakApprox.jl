module SmolyakApprox

using LinearAlgebra
using ThreadPools

include("smolyak_approx_functions.jl")

export chebyshev_gauss_lobatto,
       clenshaw_curtis_equidistant,
       smolyak_grid,
       smolyak_weights,
       smolyak_weights_threaded,
       smolyak_inverse_interpolation_matrix,
       smolyak_inverse_interpolation_matrix_threaded,
       smolyak_evaluate,
       smolyak_polynomial,
       smolyak_derivative_finite_difference,
       smolyak_derivative,
       smolyak_gradient,
       smolyak_grid_full,
       smolyak_weights_full,
       smolyak_evaluate_full,
       smolyak_derivative_full,
       smolyak_gradient_full,
       smolyak_pl_weights,
       smolyak_pl_weights_threaded,
       smolyak_pl_evaluate

end
