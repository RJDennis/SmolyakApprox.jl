# Solve the model using policy function iteration

using SmolyakApprox
using GaussQuadrature
using NLsolve

function solve_stoch_growth_model()

  function cond_future_shock(shock_state)

    cond_future_tech_shock = rho*shock_state + sqrt(2)*sigma_eps*eps_nodes

    return cond_future_tech_shock

  end

  function stoch_growth_model(x::Array{Float64,1},f::Array{Float64,1})

    # x = [capital t+1, consumption t]

    cond_future_a_shock = cond_future_shock(state[1])

    cond_expectation_consumption = 0.0
    for i = 1:number_eps_nodes
      cond_expectation_consumption += smolyak_evaluate(consumption_weights,[cond_future_a_shock[i];x[1]],multi_ind,domain)^(-sigma)*(1.0 - delta + alpha*exp(cond_future_a_shock[i])*x[1]^(alpha-1.0))*pi^(-1/2)*eps_weights[i]
    end

    f[1] = (1.0-delta)*state[2] + exp(state[1])*state[2]^alpha - x[2] - x[1]
    f[2] = (beta*cond_expectation_consumption)^(-1/sigma) - x[2]

  end

  beta      = 0.99
  sigma     = 2.00
  delta     = 0.06
  alpha     = 0.40
  rho       = 0.95
  sigma_eps = 0.01

  d  = 2       # Set the number of dimensions
  mu = [5, 5]  # Set the level, anisotropic

  k_max = 26.0
  k_min = 10.0
  k_domain = [k_max, k_min]

  a_max =  3*sqrt(sigma_eps^2/(1.0-rho^2))
  a_min = -3*sqrt(sigma_eps^2/(1.0-rho^2))
  a_domain = [a_max, a_min]

  domain = [a_domain k_domain]
  grid, multi_ind = smolyak_grid(chebyshev_gauss_lobatto,d,mu,domain)  # Construct the Smolyak grid and the multi index

  number_eps_nodes = 11
  (eps_nodes,eps_weights) = hermite(number_eps_nodes)

  capital     = zeros(size(grid,1))
  consumption = zeros(size(grid,1))
  for i = 1:size(grid,1)
    capital[i]     = 6*exp(grid[i,1])*grid[i,2]^alpha
    consumption[i] = 0.8*exp(grid[i,1])*grid[i,2]^alpha
  end
  consumption_weights = smolyak_weights(consumption,grid,multi_ind,domain)

  new_capital     = similar(capital)
  new_consumption = similar(consumption)

  state = zeros(2)

  len = Inf

  while len > 1e-5

    consumption_weights = smolyak_weights(consumption,grid,multi_ind,domain)

    for i = 1:size(grid,1)

      state = grid[i,:]
      initial = [capital[i]; consumption[i]]

      soln = nlsolve(stoch_growth_model,initial)
      new_capital[i]     = soln.zero[1]
      new_consumption[i] = soln.zero[2]

    end

    len = maximum(abs,[(new_consumption-consumption);(new_capital-capital)])
    println(len)
    capital     = copy(new_capital)
    consumption = copy(new_consumption)

  end

  return capital, consumption

end

solve_stoch_growth_model()
