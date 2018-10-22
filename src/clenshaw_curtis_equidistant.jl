function clenshaw_curtis_equidistant(n::S,domain = [1.0,-1.0]) where {S<:Integer}

  # Construct the nodes on the [-1.0,1.0] interval

  if n == 1
    nodes   = [0.0]
    weights = [2.0]
  else
    nodes    = zeros(n)
    nodes[1] = -1.0
    nodes[n] = 1.0

    weights    = zeros(n)
    weights[1] = 2.0/n
    weights[n] = 2.0/n

    for i = 2:(n-1)
      nodes[i]   = 2*(i-1)/(n-1)-1.0
      weights[i] = 2.0/n
    end

    if isodd(n)
      nodes[round(Int,(n+1)/2)] = 0.0
    end
  end

  # Scale the nodes to the desired domain

  nodes = scale_nodes(nodes,domain)

  return nodes, weights

end
