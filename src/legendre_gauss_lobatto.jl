function legendre_gauss_lobatto(n::S,domain = [1.0,-1.0]) where S <: Integer


  # Construct the nodes on the [-1.0,1.0] interval

  p = zeros(n,n)

  if n == 1

    nodes   = [0.0]
    weights = [2.0]

  else

    nodes    = zeros(n)
    nodes[1] = 1.0
    nodes[n] = -1.0

    for i = 2:(n-1)
      nodes[i] = cos(pi*(i-1)/(n-1))
    end

    if isodd(n)
      nodes[round(Int,(n+1)/2)] = 0.0
    end

    p[:,1] = 1.0

    len = Inf
    while len > eps(1.0)

      nodes_old = copy(nodes)
      p[:,2] = nodes

      for i = 2:(n-1)
        p[:,i+1] = ((2*i-1)*nodes.*p[:,i]-(i-1)*p[:,i-1])/i
      end

      nodes = nodes_old - (nodes.*p[:,n]-p[:,n-1])./(n*p[:,n])
      len = maximum(abs,nodes-nodes_old)

    end

  end

  # Scale the nodes to the desired domain

  nodes = domain[2] .+ (1.0 .+ nodes)*(domain[1]-domain[2])/2

  # Compute the weights

  weights = 2.0./((n-1)*n*p[:,n].^2)

  return nodes, weights

end
