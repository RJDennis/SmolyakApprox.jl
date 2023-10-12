module SmolyakApprox

import ChebyshevApprox: chebyshev_extrema,
                        normalize_node,
                        chebyshev_polynomial,
                        chebyshev_polynomial_deriv,
                        chebyshev_polynomial_sec_deriv

using ThreadPools
using LinearAlgebra

include("smolyak_approx_functions.jl")

export SApproxPlan

export chebyshev_gauss_lobatto,
       clenshaw_curtis_equidistant,
       smolyak_grid,
       smolyak_plan,
       smolyak_weights,
       smolyak_weights_threaded,
       smolyak_inverse_interpolation_matrix,
       smolyak_inverse_interpolation_matrix_threaded,
       smolyak_pl_weights,
       smolyak_pl_weights_threaded,
       smolyak_polynomial,
       smolyak_evaluate,
       smolyak_pl_evaluate,
       smolyak_interp,
       smolyak_interp_threaded,
       smolyak_derivative,
       smolyak_gradient,
       smolyak_gradient_threaded,
       smolyak_hessian,
       smolyak_hessian_threaded,
       smolyak_integrate,
       smolyak_clenshaw_curtis,
       smolyak_grid_full,
       smolyak_weights_full,
       smolyak_evaluate_full,
       smolyak_derivative_full,
       smolyak_gradient_full

end
