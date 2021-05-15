# This code presents an example to illustrate how SmolyakApprox can be used

using SmolyakApprox

function test_smolyak_approx()

  d  = 5  # Set the number of dimensions
  mu = 3  # Set the level of approximation

  g2, m2 =  smolyak_grid(clenshaw_curtis_equidistant,d,mu)
  grid, multi_ind = smolyak_grid(chebyshev_gauss_lobatto,d,mu)  # Construct the Smolyak grid and the multi index

  # An arbitrary test function

  function test(grid)

    y_value = (grid[:,1].+1).^0.1.*exp.(grid[:,2]).*log.(grid[:,3].+2).^0.2.*(grid[:,4].+2).^0.8.*(grid[:,5].+7).^0.1

    return y_value

  end

  y = test(grid)  # Evaluate the test function on the Smolyak grid
  y_pl = test(g2)

  point = [0.75, 0.45, 0.82, -0.15, -0.95]  # Choose a point to evaluate the approximated function

  # One way of computing the weights and evaluating the approximated function

  weights = smolyak_weights(y,grid,multi_ind)        # Compute the Smolyak weights
  y_hat = smolyak_evaluate(weights,point,multi_ind)  # Evaluate the approximated function

  #= A second way of computing the weights and evaluating the approximated function that
  computes the interpolation matrix just once. =#

  interp_mat = smolyak_inverse_interpolation_matrix(grid,multi_ind)  # Compute the interpolation matrix
  w = smolyak_weights(y,interp_mat)             # Compute the Smolyak weights
  y_hatt = smolyak_evaluate(w,point,multi_ind)  # Evaluate the approximated function

  # Piecewise linear

  w_pl = smolyak_pl_weights(y_pl,g2,m2)
  y_pl_hat = smolyak_pl_evaluate(w_pl,point,g2,m2)

  # Evaluate the exact function at point

  y_actual = test(point')

  # Now consider the ansiotropic case

  mu = [3, 2, 2, 2, 3]
  grid, multi_ind = smolyak_grid(chebyshev_gauss_lobatto,d,mu)  # Construct the Smolyak grid and the multi index
  y = test(grid)
  weights = smolyak_weights(y,grid,multi_ind)                   # Compute the Smolyak weights
  y_hat_ansio = smolyak_evaluate(weights,point,multi_ind)       # Evaluate the approximated function
  weights_th = smolyak_weights_threaded(y,grid,multi_ind)                   # Compute the Smolyak weights

  g3, m3 =  smolyak_grid(clenshaw_curtis_equidistant,d,mu)
  y_pl = test(g3)
  w_pl_ansio = smolyak_pl_weights(y_pl,g3,m3)
  y_pl_hat_ansio = smolyak_pl_evaluate(w_pl,point,g3,m3)
  w_pl_ansio_th = smolyak_pl_weights_threaded(y_pl,g3,m3)

  # Now test the full grid results

  mu = 3
  grid_full, multi_ind_full = smolyak_grid_full(chebyshev_gauss_lobatto,d,mu) # Construct the Smolyak grid and the multi index
  y_full = test(grid_full)
  weights_full = smolyak_weights_full(y_full,grid_full,multi_ind_full)        # Compute the Smolyak weights
  y_hat_full = smolyak_evaluate_full(weights_full,point,multi_ind_full)       # Evaluate the approximated function

  return y_hat, y_hatt, y_pl_hat, y_actual, y_hat_ansio, y_pl_hat_ansio, y_hat_full

end

test_smolyak_approx()
