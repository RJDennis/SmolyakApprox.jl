function chebyshev_gauss_lobatto(n::S,domain::Array{T,1} = [1.0,-1.0]) where {S<:Integer,T<:AbstractFloat}

  # These nodes a just the Chebyshev extrema by a different name.

  # Construct the nodes on the [-1.0,1.0] interval

  if n <= 0
    error("The number of nodes must be positive.")
  end

  if n == 1
    nodes   = [0.0]
    weights = [1.0*pi]
  else
    nodes    = zeros(n)
    nodes[1] = -1.0
    nodes[n] = 1.0

    weights    = zeros(n)
    weights[1] = pi/(2*(n-1))
    weights[n] = pi/(2*(n-1))

    for i = 2:div(n,2)
      x = cos(pi*(i-1)/(n-1))
      nodes[i]       = -x
      nodes[end-i+1] = x
      weights[i]       = pi/(n-1)
      weights[end-i+1] = pi/(n-1)
    end

    if isodd(n)
      nodes[round(Int,(n+1)/2)]   = 0.0
      weights[round(Int,(n+1)/2)] = pi/(n-1)
    end
  end

  # Scale the nodes to the desired domain

  nodes = scale_nodes(nodes,domain)

  return nodes, weights

end