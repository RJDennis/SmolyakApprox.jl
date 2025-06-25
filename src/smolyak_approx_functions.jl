import ChebyshevApprox: chebyshev_extrema,
                        normalize_node,
                        chebyshev_polynomial,
                        chebyshev_polynomial_deriv,
                        chebyshev_polynomial_sec_deriv

abstract type SApproximationPlan end

"""
SApproxPlan is an immutable struct that contains the information used to approximate a multi-dimensional function.
SApproxPlan has four fileds: A node type that describes how the approximating points are generated, the approximation
grid, the multi-index associated with the Smolyak polynomial, and the approximation domain.
""" 
struct SApproxPlan{S<:Integer,T<:AbstractFloat} <: SApproximationPlan

  node_type::Symbol
  grid::Union{Array{T,1},Array{T,2}}
  multi_index::Union{Array{S,1},Array{S,2}}
  domain::Union{Array{T,1},Array{T,2}}

end

"""
Computes ```N``` Chebyshev-Gauss-Lobatto (Chebyshev extrema) and scales the points to the interval given in
```domain``` (defaults to [1.0,-1.0]).  Returns a vector containing the ```N``` points located over the 
specified domain.

Chebyshev-Gauss-Lobatto points are used in the case where the Smolyak approximation is based on Chebyshev
polynomials.

Signatures
==========

points = chebyshev_gauss_lobatto(N,domain)
points = chebyshev_extrema(N,domain)

Examples
========
```
julia> points = chebyshev_gauss_lobatto(4)
[-1.0
 -0.5000000000000001
  0.5000000000000001
  1.0]

julia> points = chebyshev_gauss_lobatto(4,[3.0,-1.0])
[-1.0
 -2.220446049250313e-16
  2.0
  3.0]
```
"""
const chebyshev_gauss_lobatto = chebyshev_extrema

"""
Computes ```N``` uniformly spaced points on [1.0,-1.0] and scales the points to the interval given in ```domain```
(defaults to [1.0,-1.0]).  Returns a vector containing the ```N``` points located over the specified domain.

clenshaw-Curtis equidistant points are used in the case where the Smolyak approximation is based on piecewise
linear basis functions.

Signature
=========

points = clenshaw_curtis_equidistant(N,domain)

Examples
========
```
julia> points = clenshaw_curtis_equidistant(4)
[-1.0
 -0.33333333333333337
  0.3333333333333335
  1.0]

julia> points = clenshaw_curtis_equidistant(4,[3.0,-1.0])
[-1.0
  0.33333333333333326
  1.666666666666667
  3.0]
```
"""
function clenshaw_curtis_equidistant(n::S,domain = [1.0,-1.0]) where {S<:Integer}

  # Construct the nodes on the [-1.0,1.0] interval

  if n == 1
    nodes   = [0.0] # allocates
  else
    nodes    = zeros(n) # allocates
    nodes[1] = -1.0
    nodes[n] =  1.0
    for i = 2:div(n,2)
      nodes[i]       =  2*(i-1)/(n-1) - 1.0
      nodes[end-i+1] = -2*(i-1)/(n-1) + 1.0
    end
  end

  # Scale the nodes to the desired domain

  scale_nodes!(nodes,domain)

  return nodes

end

"""
Computes ```N``` uniformly spaced points of floating point type ```T``` on [1.0,-1.0] and scales the points to the
interval given in ```domain``` (defaults to [1.0,-1.0]).  Returns a vector containing the ```N``` points located 
over the specified domain.

Signature
=========

points = clenshaw_curtis_equidistant(T,N,domain)

Examples
========
```
julia> using DoubleFloats
julia> points = clenshaw_curtis_equidistant(Double64,4)
[-1.0
 -0.33333333333333337
  0.33333333333333337
  1.0]

julia> points = clenshaw_curtis_equidistant(Double64,4,[3.0,-1.0])
[-1.0
  0.33333333333333326
  1.6666666666666667
  3.0]

  julia> points = clenshaw_curtis_equidistant(BigFloat,4,[3.0,-1.0])
[-1.0
  0.3333333333333332593184650249895639717578887939453125
  1.6666666666666667406815349750104360282421112060546875
  3.0]
```
"""
function clenshaw_curtis_equidistant(T::DataType,n::S,domain = [1.0,-1.0]) where {S<:Integer}

  # Construct the nodes on the [-1.0,1.0] interval

  if n == 1
    nodes   = [zero(T)] # allocates
  else
    nodes    = zeros(T,n) # allocates
    nodes[1] = -one(T)
    nodes[n] =  one(T)
    for i = 2:div(n,2)
      nodes[i]       =  2*(i-1)/(n-1) - one(T)
      nodes[end-i+1] = -2*(i-1)/(n-1) + one(T)
    end
  end

  # Scale the nodes to the desired domain

  scale_nodes!(nodes,domain)

  return nodes

end

"""
Creates the multi-index for a ```d```-variable ansiotropic grid with layers given by the vector ```mu```.  Returns
a matrix orf integers.

Signature
=========

m_index = generate_multi_index(d,mu)

Example
=======
```
julia> m_index = generate_multi_index(2,[2,1])
[1  1
 2  1
 3  1
 1  2
 2  2]
"""
function generate_multi_index(d::S,mu::Array{S,1}) where {S<:Integer}

  nt = num_terms(mu,d)
  multi_index = Array{S,2}(undef,nt,d) # allocates
  multi_index[1,:] = ones(S,1,d) # allocates

  max_mu = maximum(mu)

  w = Tuple(max_mu+1 for _ in 1:d) # allocates
  candidate_indexes = Tuple.(CartesianIndices(w)) # allocates
  pos = 1
  @inbounds for i = 2:(max_mu+1)^d
    if sum(candidate_indexes[i]) <= d+max_mu && sum(candidate_indexes[i] .<= mu .+ 1) == d
      pos += 1
      if pos > nt # handles the case where nt is under-estimated
        multi_index = [multi_index; collect(candidate_indexes[i])'] # allocates
      else
       multi_index[pos,:] .= candidate_indexes[i] # allocates
      end
    end
  end

  if pos < nt # handles case where nt is over-estimated
    multi_index = multi_index[1:pos,:]
  end

  return multi_index

end

"""
Creates the multi-index for a ```d```-variable isotropic grid with layers given by the integer ```mu```. Returns
a matrix of integers.

Signature
=========

m_index = generate_multi_index(d,mu)

Example
=======
```
julia> m_index = generate_multi_index(2,2])
[1  1
 2  1
 3  1
 1  2
 2  2
 1  3]
"""
function generate_multi_index(d::S,mu::S) where {S<:Integer}

  if d == 1
    multi_index = [i for i in 1:mu+1]
    return multi_index
  else
    multi_index_base = generate_multi_index(d-1,mu)
    N = size(multi_index_base,1)
    multi_index = zeros(S,N*(mu+1),d)
    pos = 0
    @inbounds @views for j = 1:N
      for i = 1:mu+1
        if sum(multi_index_base[j,:]) + i <= d+mu
          pos += 1
          multi_index[pos,2:d] .= multi_index_base[j,:]
          multi_index[pos,1] = i
        end
      end
    end
    return multi_index[1:pos,:]
  end
end    
  
# Computes the number of terms in the multi-index for the isotropic case (it also computes the number of terms in
# a complete polynominal based on the order and the number of dimensions).

function num_terms(order::S,d::S) where {S <: Integer} # Internal function, not exported

  if d == 1
    return order+1
  else
    return div(num_terms(order,d-1)*(order+d),d)
  end

end

# Poorly approximates the number of terms in the multi-index for the ansiotropic case.

function num_terms(order::Array{S,1},d::S) where {S<:Integer} # Internal function, not exported

  max_mu = maximum(order)
  nt = num_terms(max_mu,d) # Deliberate over-estimate of the number of terms
    
  return nt
  
end

m_i(x::S) where {S <: Integer} = (x == 1 ? 1 : 2^(x-1) + 1) # Internal function, not exported
m_i(x::Array{S,N}) where {S <: Integer,N} = m_i.(x) # Internal function not exported

function combine_nodes(nodes1::Union{Array{R,1},Array{R,2}},nodes2::Array{R,1}) where {R<:Number} # Internal function, not exported
  
  # nodes1 can be a 1d or 2d array; nodes2 is a 1d array

  n1 = size(nodes1,1)
  n2 = size(nodes1,2)
  n3 = length(nodes2)

  combined_nodes = Array{R,2}(undef,n1*n3,n2+1)

  @inbounds for i = 1:n3
    combined_nodes[(i-1)*n1+1:i*n1,1:n2] = nodes1
  end
  @inbounds for i = 1:n1
    @inbounds for j = 1:n3
      combined_nodes[(j-1)*n1+i,n2+1] = nodes2[j]
    end
  end

  return combined_nodes

end

function scale_nodes!(nodes::Array{R,1},domain::Array{T,1}) where {T<:AbstractFloat,R<:Number} # Internal function, not exported

  @inbounds for i in eachindex(nodes)
    nodes[i] = domain[2] + (1.0+nodes[i])*(domain[1]-domain[2])*0.5
  end

end

function scale_nodes!(nodes::Array{R,2},domain::Array{T,2}) where {T<:AbstractFloat,R<:Number} # Internal function, not exported

  @inbounds for i in CartesianIndices(nodes)
    nodes[i] = domain[2,i[2]] + (1.0+nodes[i])*(domain[1,i[2]]-domain[2,i[2]])*0.5
  end
   
end

"""
Ensures the approximation domain has the expected form.

Signatures
==========

domain = check_domain(dom)

Example
=======

```
julia> dom = [-1.0 2.0; 2.0 -1.0]
julia> domain = check_domain(dom)
```
"""
function check_domain(dom::Array{T,2}) where {T<:AbstractFloat}

  domain = similar(dom)

  for i in axes(dom,2)
    domain[1,i] = max(dom[1,i],dom[2,i])
    domain[2,i] = min(dom[1,i],dom[2,i])
  end

  return domain

end

"""
Uses the ```node_type``` function to construt the ```d```-dimensional Smolyak grid with approximation layer
```mu``` and ```domain```.  If ```mu``` is an integer (vector of integers) the isotropic (ansiotropic) grid is
constructed.  Returns the approximation grid and the associated multi-index.  If ```domain``` is not provided,
then the approximation domain defaults to [1.0,-1.0]^d.

Signatures
==========

grid, multi_index = smolyak_grid(node_type,d,mu)
grid, multi_index = smolyak_grid(chebyshev_extrema,d,mu,domain)

Examples
========
```
julia> grid, m_index = smolyak_grid(chebyshev_extrema,2,2)
julia> grid, m_index = smolyak_grid(chebyshev_extrema,2,[2,1])
julia> grid, m_index = smolyak_grid(chebyshev_extrema,2,2,[3.0 1.5; 2.0 0.5])
julia> grid, m_index = smolyak_grid(chebyshev_extrema,2,[2,2],[3.0 1.5; 2.0 0.5])
```
"""
function smolyak_grid(node_type::Function,d::S,mu::Union{S,Array{S,1}},domain=[ones(1, d); -ones(1, d)]) where {S<:Integer}

  if size(domain,2) != d
    error("'domain' has too many or too few columns.")
  end

  T = eltype(domain)

  multi_index        = generate_multi_index(d,mu)
  unique_multi_index = sort(unique(multi_index))
  unique_node_number = m_i.(unique_multi_index)

  # Create base nodes to be used in the sparse grid

  base_nodes   = Array{Array{T,1},1}(undef,length(unique_node_number))
  for i in eachindex(unique_node_number)
    base_nodes[i] = node_type(unique_node_number[i])
  end

  # Determine the unique nodes introduced at each higher level

  unique_base_nodes = Array{Array{T,1},1}(undef,length(unique_node_number))
  unique_base_nodes[1] = base_nodes[1]
  for i = 2:length(unique_base_nodes)
    unique_base_nodes[i] = setdiff(base_nodes[i],base_nodes[i-1])
  end

  # Construct the sparse grid from the unique nodes

  grid = Array{T,2}(undef,determine_grid_size(multi_index))
  l = 1
  @inbounds for j in axes(multi_index,1)
    new_nodes = unique_base_nodes[multi_index[j,1]] # Here new_nodes is a 1d array
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,unique_base_nodes[multi_index[j,i]])  # Here new_nodes becomes a 2d array
    end
    m = size(new_nodes,1)
    grid[l:l+m-1,:] = new_nodes
    l += m
  end

  if d == 1
    grid = grid[:]
  end
  
  # Now scale the nodes to the desired domain

  scale_nodes!(grid,domain)

  return grid, multi_index

end

# Uses the multi-index to determine the grid size.

function determine_grid_size(mi::AbstractArray{S,N}) where {S<:Integer,N} # Internal function, not exported

  temp = similar(mi)

  for i in axes(mi,1)
    for j in axes(mi,2)
      if mi[i,j] == 1
        temp[i,j] = 1
      elseif mi[i,j] == 2
        temp[i,j] = 2^(mi[i,j]-1)
      else
        temp[i,j] = 2^(mi[i,j]-1)+1 - (2^(mi[i,j]-2)+1)
      end
    end
  end

  s = 0
  for i in axes(mi,1)
    t = 1
    for j in axes(mi,2)
      t *= temp[i,j]
    end
    s += t
  end

  return (s, size(mi,2))

end

"""
Constructs a Smolyak approximation plan, given the ```node_type``` function, the number of dimensions, ```d```, the
approximation layers, ```mu```, and the approximation ```domain``` (defaults to [1.0,-1.0]^d).  Returns an SApproxPlan
struct.

Signature
=========

splan = smolyal_plan(node_type,d,mu,domain)

Examples
========
```
julia> splan = smolyak_plan(chebyshev_extrema,2,2)
julia> splan = smolyak_plan(clenshaw_curtis_equidistant,2,[2,2])
julia> splan = smolyak_plan(chebyshev_extrema,2,2,[3.0 1.5; 2.0 0.5])
julia> splan = smolyak_plan(clenshaw_curtis_equidistant,2,[2,2],[3.0 1.5; 2.0 0.5])
```
"""
function smolyak_plan(node_type::Function,d::S,mu::Union{S,Array{S,1}},domain=[ones(1, d); -ones(1, d)]) where {S<:Integer}

  g, mi = smolyak_grid(node_type,d,mu,domain)

  plan = SApproxPlan(Symbol(node_type),g,mi,domain)

  return plan

end

"""
Uses Chebyshev polynomials as basis functions to compute the weights in a Smolyak polynomial approximation given the
approximation sample, ```y```, the approximation ```grid```, the ```multi_index```, and the approximation 
```domain``` (defaults to [1.0,-1.0]^d).  Returns a vector containing the weights in the Smolyak polynomial.

Signature
=========

w = smolyak_weights(y,grid,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> f(x) = sum(x.^2)
julia> y = [f(g[i,:]) for i in axes(g,1)]
julia> w = smolyak_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
```
"""
function smolyak_weights(y::Array{T,1},grid::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain=[ones(1,size(grid,2));-ones(1,size(grid,2))]) where {T<:AbstractFloat,S<:Integer}

  # Normalize grid to the [1.0,-1.0]^d interval

  grid = copy(grid)
  for i in axes(grid,2)
    grid[:,i] = normalize_node(grid[:,i],domain[:,i])
  end

  interpolation_matrix = zeros(size(grid,1),size(grid,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i.(unique_multi_index) .- 1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct a row of the interpolation matrix
  #   Iterate over the grid, doing the above for steps at each iteration, to compute all rows of the interpolation matrix

  base_polynomials        = Array{Array{T,2},1}(undef,length(unique_orders))
  unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_orders))

  @inbounds for k in axes(grid,1)

    # Construct the base polynomials

    for i in eachindex(unique_orders)
      base_polynomials[i] = chebyshev_polynomial(unique_orders[i],grid[k,:])
    end

    # Compute the unique polynomial terms from the base polynomials

    for i = length(unique_orders):-1:2
      unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    end
    unique_base_polynomials[1] = base_polynomials[1]

    # Construct a row of the interplation matrix

    l = 1
    @inbounds @views for j in axes(multi_index,1)
      new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
      end
      m = length(new_polynomials)
      interpolation_matrix[k,l:l+m-1] = new_polynomials
      l += m
    end

  end

  weights = interpolation_matrix\y

  return weights

end

"""
Uses Chebyshev polynomials as basis functions to compute using multi-threading the weights in a Smolyak polynomial
approximation given the approximation sample, ```y```, the approximation ```grid```, the ```multi_index```, and 
the approximation ```domain``` (defaults to [1.0,-1.0]^d).  Returns a vector containing the weights in the Smolyak
polynomial.

Signature
=========

w = smolyak_weights_threaded(y,grid,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> f(x) = sum(x.^2)
julia> y = [f(g[i,:]) for i in axes(g,1)]
julia> w = smolyak_weights_threaded(y,g,mi,[1.0 1.0; 0.0 0.0])
```
"""
function smolyak_weights_threaded(y::Array{T,1},grid::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain=[ones(1,size(grid,2));-ones(1,size(grid,2))]) where {T<:AbstractFloat,S<:Integer}

  # Normalize grid to the [1.0,-1.0]^d interval

  grid = copy(grid)
  for i in axes(grid,2)
    grid[:,i] = normalize_node(grid[:,i],domain[:,i])
  end

  interpolation_matrix = zeros(size(grid,1),size(grid,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i.(unique_multi_index) .- 1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct a row of the interpolation matrix
  #   Iterate over the grid, doing the above three steps at each iteration, to compute all rows of the interpolation matrix

  @inbounds @sync Threads.@threads for k in axes(grid,1)

    base_polynomials        = Array{Array{T,2},1}(undef,length(unique_orders))
    unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_orders))

    # Construct the base polynomials

    for i in eachindex(unique_orders)
      base_polynomials[i] = chebyshev_polynomial(unique_orders[i],grid[k,:])
    end

    # Compute the unique polynomial terms from the base polynomials

    for i = length(unique_orders):-1:2
      unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    end
    unique_base_polynomials[1] = base_polynomials[1]

    # Construct a row of the interplation matrix

    l = 1
    @inbounds @views for j in axes(multi_index,1)
      new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
      end
      m = length(new_polynomials)
      interpolation_matrix[k,l:l+m-1] = new_polynomials
      l += m
    end

  end

  weights = interpolation_matrix\y

  return weights

end

"""
Uses Chebyshev polynomials as basis functions to compute the weights in a Smolyak polynomial
approximation given the approximation sample, ```y```, the ```inverse interpolation matrix```.
Returns a vector containing the weights in the Smolyak polynomial.

Signature
=========

w = smolyak_weights(y,grid,inverse_interpolation_matrix)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> f(x) = sum(x.^2)
julia> y = [f(g[i,:]) for i in axes(g,1)]
julia> iim = smolyak_inverse_interpolation_matrix(g,mi)
julia> w = smolyak_weights(y,iim)
```
"""
function smolyak_weights(y::Array{T,1},inverse_interpolation_matrix::Array{T,2}) where {T<:AbstractFloat}

  weights = inverse_interpolation_matrix*y

  return weights

end

"""
Uses Chebyshev polynomials as basis functions to compute the inverse interpolation matrix for a Smolyak polynomial
approximation given the approximation grid, ```grid```, the ```multi_index```, and the approximation ```domain``` 
(defaults to [1.0,-1.0]^d).  Returns a matrix containing the inverse interoplation matrix.

Signature
=========

iim = smolyak_inverse_interpolation_matrix(grid,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> iim = smolyak_inverse_interpolation_matrix(g,mi)
```
"""
function smolyak_inverse_interpolation_matrix(grid::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain=[ones(1,size(grid,2));-ones(1,size(grid,2))]) where {T<:AbstractFloat,S<:Integer}

  # Normalize grid to the [1.0,-1.0]^d interval
  
  grid = copy(grid)
  for i in axes(grid,2)
    grid[:,i] = normalize_node(grid[:,i],domain[:,i])
  end

  interpolation_matrix = zeros(size(grid,1),size(grid,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i.(unique_multi_index) .- 1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct a row of the interpolation matrix
  #   Iterate over the grid, doing the above three steps at each iteration, to compute all rows of the interpolation matrix

  base_polynomials        = Array{Array{T,2},1}(undef,length(unique_orders))
  unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_orders))

  @inbounds for k in axes(grid,1)

    # Construct the base polynomials

    for i in eachindex(unique_orders)
      base_polynomials[i] = chebyshev_polynomial(unique_orders[i],grid[k,:])
    end

    # Compute the unique polynomial terms from the base polynomials

    for i = length(unique_orders):-1:2
      unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    end
    unique_base_polynomials[1] = base_polynomials[1]

    # Construct the first row of the interplation matrix

    l = 1
    @inbounds @views for j in axes(multi_index,1)
      new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
      end
      m = length(new_polynomials)
      interpolation_matrix[k,l:l+m-1] = new_polynomials
      l += m
    end

  end

  inverse_interpolation_matrix = inv(interpolation_matrix)

  return inverse_interpolation_matrix

end

"""
Uses Chebyshev polynomials as basis functions to compute using multi-threading the inverse interpolation matrix for
a Smolyak polynomial approximation given the approximation grid, ```grid```, the ```multi_index```, and the 
approximation ```domain``` (defaults to [1.0,-1.0]^d).  Returns a matrix containing the inverse interoplation matrix.

Signature
=========

iim = smolyak_inverse_interpolation_matrix_threaded(grid,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> iim = smolyak_inverse_interpolation_matrix_threaded(g,mi)
```
"""
function smolyak_inverse_interpolation_matrix_threaded(grid::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain=[ones(1,size(grid,2));-ones(1,size(grid,2))]) where {T<:AbstractFloat,S<:Integer}

  # Normalize grid to the [1.0,-1.0]^d interval
    
  grid = copy(grid)
  for i in axes(grid,2)
    grid[:,i] = normalize_node(grid[:,i],domain[:,i])
  end

  interpolation_matrix = zeros(size(grid,1),size(grid,1))

  unique_multi_index = sort(unique(multi_index))
  unique_orders      = m_i.(unique_multi_index) .- 1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct a row of the interpolation matrix
  #   Iterate over the grid, doing the above three steps at each iteration, to compute all rows of the interpolation matrix

  @inbounds @sync Threads.@threads for k in axes(grid,1)

    base_polynomials        = Array{Array{T,2},1}(undef,length(unique_orders))
    unique_base_polynomials = Array{Array{T,2},1}(undef,length(unique_orders))

    # Construct the base polynomials

    for i in eachindex(unique_orders)
      base_polynomials[i] = chebyshev_polynomial(unique_orders[i],grid[k,:])
    end

    # Compute the unique polynomial terms from the base polynomials

    for i = length(unique_orders):-1:2
      unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    end
    unique_base_polynomials[1] = base_polynomials[1]

    # Construct the first row of the interplation matrix

    l = 1
    @inbounds @views for j in axes(multi_index,1)
      new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
      end
      m = length(new_polynomials)
      interpolation_matrix[k,l:l+m-1] = new_polynomials
      l += m
    end

  end

  inverse_interpolation_matrix = inv(interpolation_matrix)

  return inverse_interpolation_matrix

end

"""
Computes the weights in a Smolyak polynomial approximation formed using piecewise linear basis functions given the 
approximation sample, ```y```, the approximation ```grid```, the ```multi_index```, and the approximation 
```domain``` (defaults to [1.0,-1.0]^d).  Returns a vector containing the weights in the Smolyak
polynomial.

Signature
=========

w = smolyak_pl_weights(y,grid,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid(clenshaw_curtis_equidistant,2,2,[1.0 1.0; 0.0 0.0])
julia> f(x) = sum(x.^2)
julia> y = [f(g[i,:]) for i in axes(g,1)]
julia> w = smolyak_pl_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
```
"""
function smolyak_pl_weights(y::AbstractArray{T,1},grid::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain=[ones(1,size(grid,2));-ones(1,size(grid,2))]) where {T<:AbstractFloat,S<:Integer}

  # Normalize grid to the [1.0,-1.0]^d interval
    
  grid = copy(grid)
  for i in axes(grid,2)
    grid[:,i] = normalize_node(grid[:,i],domain[:,i])
  end

  interpolation_matrix = zeros(size(grid,1),size(grid,1))

  @inbounds for l in axes(grid,1)
    k = 1
    x = grid[l,:]
    @inbounds for i in axes(multi_index,1)
      m_node_number = m_i.(multi_index[i,:])
      if prod(m_node_number) == 1
        interpolation_matrix[l,k] = one(T)
        k += 1
      else
        extra_nodes = 1
        @inbounds for j in eachindex(m_node_number)
          if m_node_number[j] > 1
            extra_nodes *= m_node_number[j] - m_i(multi_index[i,j] - 1)
          end
        end
        for h = 1:extra_nodes
          a = one(T)
          @inbounds for j in eachindex(m_node_number)
            if m_node_number[j] > 1
              if abs(x[j] - grid[k,j]) > 2/(m_node_number[j]-1)
                a *= zero(T)
              else
                a *= one(T) - ((m_node_number[j]-1)/2)*abs(x[j]-grid[k,j])
              end
            end
          end
          interpolation_matrix[l,k] = a
          k += 1
        end
      end
    end
  end

  weights = interpolation_matrix\y

  return weights

end

"""
Uses multi-threading to compute the weights in a Smolyak polynomial approximation formed using piecewise linear
basis functions given the approximation sample, ```y```, the approximation ```grid```, the ```multi_index```, 
and the approximation ```domain``` (defaults to [-1.0,1.0]^d).  Returns a vector containing the weights in the 
Smolyak polynomial.

Signature
=========

w = smolyak_pl_weights_threaded(y,grid,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid(clenshaw_curtis_equidistant,2,2,[1.0 1.0; 0.0 0.0])
julia> f(x) = sum(x.^2)
julia> y = [f(g[i,:]) for i in axes(g,1)]
julia> w = smolyak_pl_weights_threaded(y,g,mi,[1.0 1.0; 0.0 0.0])
```
"""
function smolyak_pl_weights_threaded(y::AbstractArray{T,1},grid::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain=[ones(1,size(grid,2));-ones(1,size(grid,2))]) where {T<:AbstractFloat,S<:Integer}

  # Normalize grid to the [1.0,-1.0]^d interval
      
  grid = copy(grid)
  for i in axes(grid,2)
    grid[:,i] = normalize_node(grid[:,i],domain[:,i])
  end
  
  interpolation_matrix = zeros(size(grid,1),size(grid,1))

  @inbounds @sync Threads.@threads for l in axes(grid,1)
    k = 1
    x = grid[l,:]
    @inbounds for i in axes(multi_index,1)
      m_node_number = m_i.(multi_index[i,:])
      if prod(m_node_number) == 1
        interpolation_matrix[l,k] = 1.0
        k += 1
      else
        extra_nodes = 1
        @inbounds for j in eachindex(m_node_number)
          if m_node_number[j] > 1
            extra_nodes *= m_node_number[j] - m_i(multi_index[i,j] - 1)
          end
        end
        for h = 1:extra_nodes
          a = 1.0
          @inbounds for j in eachindex(m_node_number)
            if m_node_number[j] > 1
              if abs(x[j] - grid[k,j]) > 2/(m_node_number[j]-1)
                a *= 0.0
              else
                a *= 1.0 - ((m_node_number[j]-1)/2)*abs(x[j]-grid[k,j])
              end
            end
          end
          interpolation_matrix[l,k] = a
          k += 1
        end
      end
    end
  end

  weights = interpolation_matrix\y

  return weights

end

"""
Computes a Smolyak polynomial given the ```node``` in the approximation grid, the ```multi-index```, and the approximation ```domain``` (defaults to [-1.0,1.0]^d).
Returns a vector containing the basis functions in the polynomial evaluated at ```node```.

Signature
=========

spoly = smolyak_polynomial(node,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> node = g[5,:]
julia> sploy = smolyak_polynomial(node,mi,[1.0 1.0; 0.0 0.0])
```
"""
function smolyak_polynomial(node::AbstractArray{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain=[ones(1,length(node));-ones(1,length(node))]) where {R<:Number,S<:Integer}

  # Normalize grid to the [1.0,-1.0]^d interval
      
  node = copy(node)
  for i in eachindex(node)
    node[i] = normalize_node(node[i],domain[:,i])
  end
  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1

  # Below we do the following things:

  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct the Smolyak polynomial

  # Here we construct the base polynomials

  base_polynomials = Array{Array{R,2}}(undef,length(unique_orders))
  for i in eachindex(unique_orders)
    base_polynomials[i] = chebyshev_polynomial(unique_orders[i],node)
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomials = Array{Array{R,2}}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
  end
  unique_base_polynomials[1] = base_polynomials[1]

  # Construct the first row of the interplation matrix

  n = determine_grid_size(multi_index)
  polynomial = Array{R,1}(undef,n[1])

  # Iterate over grid, doing the above three steps at each iteration

  l = 1
  @inbounds for j in axes(multi_index,1)
    new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
    for i = 2:size(multi_index,2)
      new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
    end
    m = length(new_polynomials)
    polynomial[l:l+m-1] = new_polynomials
    l += m
  end

  return polynomial

end

"""
Evaluates a Smolyak polynomial formed using Chebyshev basis functions, given the ```weights```, the ```point``` at
which to evaluate the polynomial, the ```multi_index```, and the approximation ```domain``` (defaults to [1.0,-1.0]^d).
Returns a scalar.

Signature
=========

yhat = smolyak_evaluate(weights,point,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = smolyak_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> yhat = smolyak_evaluate(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0])
0.5953026581237828
```
"""
function smolyak_evaluate(weights::Array{T,1},point::AbstractArray{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain=[ones(1,length(point));-ones(1,length(point))]) where {T<:AbstractFloat,R<:Number,S<:Integer}

  # Normalize point to the [1.0,-1.0]^d interval

  point = copy(point)
  for i in eachindex(point)
    point[i] = normalize_node(point[i],domain[:,i])
  end

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1
  
  # Below we do the following things:
  
  #   Generate the polynomial terms for each order
  #   Generate the unique polynomial terms introduced at each higher order
  #   Combine the polynomial terms to construct the Smolyak polynomial
  
  # Here we construct the base polynomials
  
  base_polynomials = Array{Array{R,2}}(undef,length(unique_orders))
  for i in eachindex(unique_orders)
    base_polynomials[i] = chebyshev_polynomial(unique_orders[i],point)
  end
  
  # Compute the unique polynomial terms from the base polynomials
  
  unique_base_polynomials = Array{Array{R,2}}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
  end
  unique_base_polynomials[1] = base_polynomials[1]
  
  # Construct the first row of the interplation matrix
  
  polynomials = Array{R,1}(undef,length(weights))
  
  # Iterate over grid, doing the above three steps at each iteration
  
  l = 1
  @inbounds for j in axes(multi_index,1)
    new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
    for i = 2:size(multi_index,2)
      new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
    end
    m = length(new_polynomials)
    polynomials[l:l+m-1] = new_polynomials
    l += m
  end
  
  estimate = zero(T)
  for i in eachindex(polynomials)
    estimate += polynomials[i]*weights[i]
  end
  
  return estimate
  
end

"""
Evaluates a Smolyak polynomial formed using Chebyshev basis functions, given the ```weights``` and a Smolyak polynomial.
Returns a scalar.

Signature
=========

yhat = smolyak_evaluate(weights,polynomial)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = smolyak_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> p = smolyak_polynomial([0.37,0.71],mi,[1.0 1.0; 0.0 0.0])
julia> yhat = smolyak_evaluate(w,p)
0.5953026581237828
```
"""
function smolyak_evaluate(weights::Array{T,1},polynomial::Array{R,1}) where {T<:AbstractFloat,R<:Number}

  estimate = weights'polynomial
  
  return estimate

end

"""
Evaluates a Smolyak polynomial formed using piecewise linear basis functions, given the ```weights```, the ```point``` at
which to evaluate the polynomial, approximation ```grid```, the ```multi_index```, and the approximation ```domain``` 
(defaults to [1.0,-1.0]^d).  Returns a scalar.

Signature
=========

yhat = smolyak_pl,evaluate(weights,point,grid,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid(clenshaw_curtis_equidistant,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = smolyak_pl_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> yhat = smolyak_pl_evaluate(w,[0.37,0.71],g,mi,[1.0 1.0; 0.0 0.0])
0.5549321821467206
```
"""
function smolyak_pl_evaluate(weights::Array{T,1},point::Array{R,1},grid::Union{Array{T,1},Array{T,2}},multi_index::Union{Array{S,1},Array{S,2}},domain=[ones(1,length(point)); -ones(1,length(point))]) where {T<:AbstractFloat,R<:Number,S<:Integer}

  grid  = copy(grid)
  point = copy(point)
  @inbounds for i in eachindex(point)
    grid[:,i] = normalize_node(grid[:,i],domain[:,i])
    point[i] = normalize_node(point[i],domain[:,i])
  end

  basis = Array{R,1}(undef,size(grid,1))
  k = 1
  @inbounds for i in axes(multi_index,1)
    m_node_number = m_i.(multi_index[i,:])
    if prod(m_node_number) == 1
      basis[k] = one(R)
      k += 1
    else
      extra_nodes = 1
      @inbounds for j in eachindex(m_node_number)
        if m_node_number[j] > 1
          extra_nodes *= m_node_number[j] - m_i(multi_index[i,j] - 1)
        end
      end
      @inbounds for h = 1:extra_nodes
        a = 1.0
        @inbounds for j in eachindex(m_node_number)
          if m_node_number[j] > 1
            if abs(point[j] - grid[k,j]) > 2/(m_node_number[j]-1)
              a *= zero(R)
            else
              a *= one(R) - ((m_node_number[j]-1)/2)*abs(point[j]-grid[k,j])
            end
          end
        end
        basis[k] = a
        k += 1
      end
    end
  end

  estimate = zero(R)
  @inbounds for i in eachindex(basis)
    estimate += basis[i]*weights[i]
  end

  return estimate

end

"""
Creates an interpolating function for a Smolyak approximation given the sampling points, ```y```, and the approximation ```plan```.
Returns an interpolating function.

Signature
=========

f = smolyak_interp(y,plan)

Examples
========
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> splan = smolyak_plan(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> f = smolyak_interp(y,splan)
julia> f([0.37,0.71])
0.5953026581237828

julia> g,mi = smolyak_grid(clenshaw_curtis_equidistant,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> splan = smolyak_plan(clenshaw_curtis_equidistant,2,2,[1.0 1.0; 0.0 0.0])
julia> f = smolyak_interp(y,splan)
julia> f([0.37,0.71])
0.5549321821467206
```
"""
function smolyak_interp(y::Array{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}

  if plan.node_type == :chebyshev_extrema || plan.node_type == :chebyshev_gauss_lobatto
    weights = smolyak_weights(y,plan.grid,plan.multi_index,plan.domain)
  elseif plan.node_type == :clenshaw_curtis_equidistant
    weights = smolyak_pl_weights(y,plan.grid,plan.multi_index,plan.domain)
  end

  function interp(x::Array{R,1}) where {R<:Number}

    if plan.node_type == :chebyshev_extrema || plan.node_type == :chebyshev_gauss_lobatto
      return smolyak_evaluate(weights,x,plan.multi_index,plan.domain)
    elseif plan.node_type == :clenshaw_curtis_equidistant
      return smolyak_pl_evaluate(weights,x,plan.grid,plan.multi_index,plan.domain)
    end

  end

  return interp

end

"""
Creates an interpolating function that uses multi-threading to construct an interpolating function given the
sampling points, ```y```, and the approximation ```plan```.  Returns an interpolating function.

Signature
=========

f = smolyak_interp_threaded(y,plan)

Examples
========
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> splan = smolyak_plan(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> f = smolyak_interp_threaded(y,splan)
julia> f([0.37,0.71])
0.5953026581237828

julia> g,mi = smolyak_grid(clenshaw_curtis_equidistant,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> splan = smolyak_plan(clenshaw_curtis_equidistant,2,2,[1.0 1.0; 0.0 0.0])
julia> f = smolyak_interp_threaded(y,splan)
julia> f([0.37,0.71])
0.5549321821467206
```
"""
function smolyak_interp_threaded(y::Array{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}

  if plan.node_type == :chebyshev_extrema || plan.node_type == :chebyshev_gauss_lobatto
    weights = smolyak_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  elseif plan.node_type ==:clenshaw_curtis_equidistant
    weights = smolyak_pl_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  end

  function interp(x::Array{R,1}) where {R<:Number}

    if plan.node_type == :chebyshev_extrema || plan.node_type == :chebyshev_gauss_lobatto
      return smolyak_evaluate(weights,x,plan.multi_index,plan.domain)
    elseif plan.node_type == :clenshaw_curtis_equidistant
      return smolyak_pl_evaluate(weights,x,plan.grid,plan.multi_index,plan.domain)
    end

  end

  return interp

end

function _smolyak_derivative(weights::Array{T,1},point::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer} # Internal function, not exported

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1

  # Here we construct the base polynomials

  base_polynomials            = Array{Array{R,2},1}(undef,length(unique_orders))
  base_polynomial_derivatives = Array{Array{R,2},1}(undef,length(unique_orders))
  for i in eachindex(unique_orders)
    base_polynomials[i]            = chebyshev_polynomial(unique_orders[i],point)
    base_polynomial_derivatives[i] = chebyshev_polynomial_deriv(unique_orders[i],point)
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomials = Array{Array{R,2},1}(undef,length(unique_orders))
  unique_base_polynomial_derivatives = Array{Array{R,2},1}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    unique_base_polynomial_derivatives[i] = base_polynomial_derivatives[i][:,size(base_polynomial_derivatives[i-1],2)+1:end]
  end
  unique_base_polynomials[1] = base_polynomials[1]
  unique_base_polynomial_derivatives[1] = base_polynomial_derivatives[1]

  # Construct the first row of the interplation matrix

  polynomials = Array{R,1}(undef,length(weights))

  # Iterate over nodes, doing the above three steps at each iteration

  l = 1
  @inbounds for j in axes(multi_index,1)
    if pos == 1
      new_polynomials = unique_base_polynomial_derivatives[multi_index[j,1]][1,:]
    else
      new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
    end
    for i = 2:size(multi_index,2)
      if pos == i
        new_polynomials = kron(new_polynomials,unique_base_polynomial_derivatives[multi_index[j,i]][i,:])
      else
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
      end
    end
    m = length(new_polynomials)
    polynomials[l:l+m-1] = new_polynomials
    l += m
  end

  evaluated_derivative = zero(T)

  for i in eachindex(polynomials)
    evaluated_derivative += polynomials[i]*weights[i]
  end

  return evaluated_derivative

end

"""
Computes the partial derivative of a Smolyak polynomial formed using Chebyshev basis functions, given the 
```weights```, the ```point``` at which to evaluate the polynomial, the ```multi_index```, the approximation 
```domain```, and the index of the variable to differentiate with respect to, ```pos```.
Returns a scalar.

Signature
=========

deriv = smolyak_derivative(weights,point,multi_index,domain,pos)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = smolyak_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> deriv1 = smolyak_derivative(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0],1)
0.40368169538532706
julia> deriv2 = smolyak_derivative(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0],2)
0.5250627466157657
```
"""
function smolyak_derivative(weights::Array{T,1},point::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}

  point = copy(point)

  for i in eachindex(point)
    point[i] = normalize_node(point[i],domain[:,i])
  end

  evaluated_derivative = _smolyak_derivative(weights,point,multi_index,pos)

  return evaluated_derivative*(2.0/(domain[1,pos]-domain[2,pos]))

end

"""
Computes the gradient a Smolyak polynomial formed using Chebyshev basis functions, given the ```weights```, the
```point``` at which to evaluate the polynomial, the ```multi_index```, and the approximation ```domain```.  
Returns a one-row matrix.

Signature
=========

grad = smolyak_gradient(weights,point,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = smolyak_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> grad = smolyak_gradient(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0])
[0.403682  0.525063]
```
"""
function smolyak_gradient(weights::Array{T,1},point::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain=[ones(1,length(point));-ones(1,length(point))]) where {T<:AbstractFloat,R<:Number,S<:Integer}

  d = length(point)
  gradient = Array{R,2}(undef,1,d)

  for i = 1:d
    gradient[i] = smolyak_derivative(weights,point,multi_index,domain,i)
  end

  return gradient

end

"""
Creates an interpolating function to compute the gradient of a Smolyak polynomial given the sampling points, ```y```,
and the approximation ```plan```.  Returns an interpolating function.

Signature
=========

grad = smolyak_gradient(y,plan)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> splan = smolyak_plan(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> grad = smolyak_gradient(y,splan)
julia> grad([0.37,0.71])
[0.403682  0.525063]
```
"""
function smolyak_gradient(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}

  if plan.node_type == :clenshaw_curtis_equidistant
    error("Not implemented for clenshaw_curtis_equidistant nodes")
  end

  weights = smolyak_weights(y,plan.grid,plan.multi_index,plan.domain)
  
  function smolyak_grad(x::Array{R,1}) where {R<:Number}
  
    return smolyak_gradient(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return smolyak_grad
  
end

"""
Creates an interpolating function that uses multi-treading to compute the gradient of a Smolyak polynomial given the
sampling points, ```y```, and the approximation ```plan```.  Returns an interpolating function.

Signature
=========

grad = smolyak_gradient_threaded(y,plan)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> splan = smolyak_plan(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> grad = smolyak_gradient_threaded(y,splan)
julia> grad([0.37,0.71])
[0.403682  0.525063]
```
"""
function smolyak_gradient_threaded(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}
  
  if plan.node_type == :clenshaw_curtis_equidistant
    error("Not implemented for clenshaw_curtis_equidistant nodes")
  end

  weights = smolyak_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  
  function smolyak_grad(x::Array{R,1}) where {R<:Number}
  
    return smolyak_gradient(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return smolyak_grad
  
end

"""
Computes the hessian a Smolyak polynomial formed using Chebyshev basis functions, given the ```weights```, the
```point``` at which to evaluate the polynomial, the ```multi_index```, and the approximation ```domain```.  
Returns a matrix.

Signature
=========

hess = smolyak_hessian(weights,point,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = smolyak_weights(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> hess = smolyak_hessian(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0])
[ -2.53475  1.06753
   1.06753  0.199234]
```
"""
function smolyak_hessian(weights::Array{T,1},point::Array{R,1},multi_index::Union{Array{S,1},Array{S,2}},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  point = copy(point)

  d = length(point)
  for i = 1:d
    point[i] = normalize_node(point[i],domain[:,i])
  end

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1

  hess = Array{T,2}(undef,d,d)

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1

  # Here we construct the base polynomials

  base_polynomials                = Array{Array{R,2},1}(undef,length(unique_orders))
  base_polynomial_derivatives     = Array{Array{R,2},1}(undef,length(unique_orders))
  base_polynomial_sec_derivatives = Array{Array{R,2},1}(undef,length(unique_orders))
  for i in eachindex(unique_orders)
    base_polynomials[i]                = chebyshev_polynomial(unique_orders[i],point)
    base_polynomial_derivatives[i]     = chebyshev_polynomial_deriv(unique_orders[i],point)
    base_polynomial_sec_derivatives[i] = chebyshev_polynomial_sec_deriv(unique_orders[i],point)
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomials                = Array{Array{R,2},1}(undef,length(unique_orders))
  unique_base_polynomial_derivatives     = Array{Array{R,2},1}(undef,length(unique_orders))
  unique_base_polynomial_sec_derivatives = Array{Array{R,2},1}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomials[i]                = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
    unique_base_polynomial_derivatives[i]     = base_polynomial_derivatives[i][:,size(base_polynomial_derivatives[i-1],2)+1:end]
    unique_base_polynomial_sec_derivatives[i] = base_polynomial_sec_derivatives[i][:,size(base_polynomial_sec_derivatives[i-1],2)+1:end]
  end
  unique_base_polynomials[1]                = base_polynomials[1]
  unique_base_polynomial_derivatives[1]     = base_polynomial_derivatives[1]
  unique_base_polynomial_sec_derivatives[1] = base_polynomial_sec_derivatives[1]

  # Construct the first row of the interplation matrix

  polynomials = Array{R,1}(undef,length(weights))

  # Iterate over nodes, doing the above three steps at each iteration

  @inbounds for c in CartesianIndices(hess)
    l = 1
    @inbounds for j in axes(multi_index,1)
      if 1 == c[1] == c[2]
        new_polynomials = unique_base_polynomial_sec_derivatives[multi_index[j,1]][1,:]
      elseif 1 == c[1] || 1 == c[2]
        new_polynomials = unique_base_polynomial_derivatives[multi_index[j,1]][1,:]
      else
        new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
      end
      for i = 2:size(multi_index,2)
        if i == c[1] == c[2]
          new_polynomials = kron(new_polynomials,unique_base_polynomial_sec_derivatives[multi_index[j,i]][i,:])
        elseif i == c[1] || i == c[2]
          new_polynomials = kron(new_polynomials,unique_base_polynomial_derivatives[multi_index[j,i]][i,:])
        else
          new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
        end
      end
      m = length(new_polynomials)
      polynomials[l:l+m-1] = new_polynomials
      l += m
    end

    evaluated_derivative = zero(T)

    for i in eachindex(polynomials)
      evaluated_derivative += polynomials[i]*weights[i]
    end

    hess[c] = evaluated_derivative*(2.0/(domain[1,c[1]]-domain[2,c[1]]))*(2.0/(domain[1,c[2]]-domain[2,c[2]]))

  end

  return hess

end

"""
Creates an interpolating function to compute the hessian of a Smolyak polynomial given the sampling points, ```y```,
and the approximation ```plan```.  Returns an interpolating function.

Signature
=========

hess = smolyak_hessian(y,plan)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> splan = smolyak_plan(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> hess = smolyak_hessian(y,splan)
julia> hess([0.37,0.71])
[-2.53475  1.06753
  1.06753  0.199234]
```
"""
function smolyak_hessian(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}

  if plan.node_type == :clenshaw_curtis_equidistant
    error("Not implemented for clenshaw_curtis_equidistant nodes")
  end

  weights = smolyak_weights(y,plan.grid,plan.multi_index,plan.domain)
  
  function smolyak_hess(x::Array{R,1}) where {R<:Number}
  
    return smolyak_hessian(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return smolyak_hess
  
end
"""
Creates an interpolating function that uses multi-threading to compute the hessian of a Smolyak polynomial given the
sampling points, ```y```, and the approximation ```plan```.  Returns an interpolating function.

Signature
=========

hess = smolyak_hessian_threaded(y,plan)

Example
=======
```
julia> g,mi = smolyak_grid(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> splan = smolyak_plan(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> hess = smolyak_hessian_threaded(y,splan)
julia> hess([0.37,0.71])
[-2.53475  1.06753
  1.06753  0.199234]
```
"""
function smolyak_hessian_threaded(y::AbstractArray{T,1},plan::P) where {T<:AbstractFloat,P<:SApproximationPlan}
  
  if plan.node_type == :clenshaw_curtis_equidistant
    error("Not implemented for clenshaw_curtis_equidistant nodes")
  end

  weights = smolyak_weights_threaded(y,plan.grid,plan.multi_index,plan.domain)
  
  function smolyak_hess(x::Array{R,1}) where {R<:Number}
  
    return smolyak_hessian(weights,x,plan.multi_index,plan.domain)
  
  end
  
  return smolyak_hess
  
end

function integrate_cheb_polys(order::S) where {S <: Integer} # Internal function, not exported

  # Integrates Chebyshev polynomials over the domain [-1,1]

  p = zeros(order+1)
    
  for i in 1:order+1
    if i == 2
      p[i] = 0.0
    else
      p[i] = ((-1)^(i-1)+1)/(1-(i-1)^2)
    end
  end

  return p

end

# Docstrings are done up to here.

"""
Numerically integrates a function, ```f```, over all dimensions by approximating the function with a Smolyak
polynomial according to the approximation ```plan```, using either :clenshaw_curtis or gauss_Chebyshev_quad 
as ```method```.  Returns a scalar.

Signature
=========

integral = smolyak_integrate(f,plan,method)

Example
=======
```
julia> f(x) = sum(x.^2)
julia> splan = smolyak_plan(chebyshev_extrema,4,8,[ones(1,4); zeros(1,4)])
julia> integral = smolyak_integrate(f,splan,:clenshaw_curtis)
1.3333333333333328
julia> integral = smolyak_integrate(f,splan,:gauss_chebyshev_quad)
1.3318637929187602
```
"""
function smolyak_integrate(f::Function,plan::SApproxPlan,method::Symbol)

  if method == :clenshaw_curtis
    integral = smolyak_clenshaw_curtis(f,plan)
  elseif method == :gauss_chebyshev_quad
    integral = smolyak_gauss_chebyshev_quad(f,plan)
  else
    error("Integration not implemented for that method")
  end

  return integral

end

"""
Uses the Clenshaw-Curtis method to numerically integrates a function, ```f```, over all dimensions by approximating
the function with a Smolyak polynomial according to the approximation ```plan```.  Returns a scalar.

Signature
=========

integral = smolyak_clenshaw_curtis(f,plan)

Example
=======
```
julia> f(x) = sum(x.^2)
julia> splan = smolyak_plan(chebyshev_extrema,4,8,[ones(1,4); zeros(1,4)])
julia> integral = smolyak_clenshaw_curtis(f,splan)
1.3333333333333328
```
"""
function smolyak_clenshaw_curtis(f::Function,plan::SApproxPlan)

  grid        = plan.grid
  multi_index = plan.multi_index
  domain      = plan.domain

  y = zeros(size(grid,1))
  for i in eachindex(y)
    y[i] = f(grid[i,:])
  end

  weights = smolyak_weights(y,grid,multi_index,domain)

  # Uses Clenshaw-Curtis to integrate over all dimensions

  unique_multi_index = sort(unique(multi_index))
  unique_orders = m_i.(unique_multi_index) .- 1

  # Here we construct the base polynomials

  T = eltype(grid)

  base_polynomial_integrals = Array{Array{T,1},1}(undef,length(unique_orders))
  for i in eachindex(unique_orders)
    base_polynomial_integrals[i] = integrate_cheb_polys(unique_orders[i])
  end

  # Compute the unique polynomial terms from the base polynomials

  unique_base_polynomial_integrals = Array{Array{T,1},1}(undef,length(unique_orders))
  for i = length(unique_orders):-1:2
    unique_base_polynomial_integrals[i] = base_polynomial_integrals[i][length(base_polynomial_integrals[i-1])+1:end]
  end
  unique_base_polynomial_integrals[1] = base_polynomial_integrals[1]

  # Construct the first row of the interplation matrix

  polynomials = Array{T,1}(undef,length(weights))

  # Iterate over nodes, doing the above three steps at each iteration

  l = 1
  @inbounds for j in axes(multi_index,1)
    new_polynomials = unique_base_polynomial_integrals[multi_index[j,1]][:]
    for i = 2:size(multi_index,2)
      new_polynomials = kron(new_polynomials,unique_base_polynomial_integrals[multi_index[j,i]][:])
    end
    m = length(new_polynomials)
    polynomials[l:l+m-1] = new_polynomials
    l += m
  end

  evaluated_integral = zero(T)

  for i in eachindex(polynomials)
    evaluated_integral += polynomials[i]*weights[i]
  end

  scale_factor = (domain[1,1]-domain[2,1])/2
  for i in 2:size(multi_index,2)
    scale_factor = scale_factor*(domain[1,i]-domain[2,i])/2
  end

  return evaluated_integral*scale_factor

end
"""
Uses the Clenshaw-Curtis method to numerically integrates a function, ```f```, over all dimensions except ```pos```
by approximating the function with a Smolyak polynomial according to the approximation ```plan```.  Returns a function.

Signature
=========

integral = smolyak_clenshaw_curtis(f,plan,pos)

Example
=======
```
julia> f(x) = sum(x.^2)
julia> splan = smolyak_plan(chebyshev_extrema,4,8,[ones(1,4); zeros(1,4)])
julia> integral = smolyak_clenshaw_curtis(f,splan,4)
julia> integral(0.5)
1.2499999999999996
```
"""
function smolyak_clenshaw_curtis(f::Function,plan::SApproxPlan,pos::S) where {S<:Integer}

  # Uses Clenshaw-Curtis to integrate over all dimensions except for pos

  grid        = plan.grid
  multi_index = plan.multi_index
  domain      = plan.domain

  y = zeros(size(grid,1))
  for i in eachindex(y)
    y[i] = f(grid[i,:])
  end

  weights = smolyak_weights(y,grid,multi_index,domain)

  function smolyak_int(point::R) where {R <: Number}
    
    point = normalize_node(point,domain[:,pos])

    unique_multi_index = sort(unique(multi_index))
    unique_orders = m_i.(unique_multi_index) .- 1

    # Here we construct the base polynomials

    T = eltype(grid)

    base_polynomials          = Array{Array{T,1},1}(undef,length(unique_orders))
    base_polynomial_integrals = Array{Array{T,1},1}(undef,length(unique_orders))
    for i in eachindex(unique_orders)
      base_polynomials[i]          = chebyshev_polynomial(unique_orders[i],point)[:]
      base_polynomial_integrals[i] = integrate_cheb_polys(unique_orders[i])
    end

    # Compute the unique polynomial terms from the base polynomials

    unique_base_polynomials          = Array{Array{T,1},1}(undef,length(unique_orders))
    unique_base_polynomial_integrals = Array{Array{T,1},1}(undef,length(unique_orders))
    for i = length(unique_orders):-1:2
      unique_base_polynomials[i]          = base_polynomials[i][length(base_polynomials[i-1])+1:end]
      unique_base_polynomial_integrals[i] = base_polynomial_integrals[i][length(base_polynomial_integrals[i-1])+1:end]
    end
    unique_base_polynomials[1]          = base_polynomials[1]
    unique_base_polynomial_integrals[1] = base_polynomial_integrals[1]

    # Construct the first row of the interplation matrix

    polynomials = Array{T,1}(undef,length(weights))

    # Iterate over nodes, doing the above three steps at each iteration

    l = 1
    @inbounds for j in axes(multi_index,1)
      if pos == 1
        new_polynomials = unique_base_polynomials[multi_index[j,1]][:]
      else
        new_polynomials = unique_base_polynomial_integrals[multi_index[j,1]][:]
      end
      for i = 2:size(multi_index,2)
        if pos == i
          new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][:])
        else
          new_polynomials = kron(new_polynomials,unique_base_polynomial_integrals[multi_index[j,i]][:])
        end
      end
      m = length(new_polynomials)
      polynomials[l:l+m-1] = new_polynomials
      l += m
    end

    evaluated_integral = zero(T)

    for i in eachindex(polynomials)
      evaluated_integral += polynomials[i]*weights[i]
    end

    scale_factor = 1.0
    for i in 1:size(multi_index,2)
      if pos != i
        scale_factor = scale_factor*(domain[1,i]-domain[2,i])/2
      end
    end

    return evaluated_integral*scale_factor

  end

  return smolyak_int

end

"""
Uses Gauss-Chebyshev quadrature to numerically integrates a function, ```f```, over all dimensions by approximating
the function with a Smolyak polynomial according to the approximation ```plan```.  Returns a scalar.

Signature
=========

integral = smolyak_gauss_chebyshev_quad(f,plan)

Example
=======
```
julia> f(x) = sum(x.^2)
julia> splan = smolyak_plan(chebyshev_extrema,4,8,[ones(1,4); zeros(1,4)])
julia> integral = smolyak_gauss_chebyshev_quad(f,splan)
1.3318637929187602
```
"""
function smolyak_gauss_chebyshev_quad(f::Function,plan::SApproxPlan)

  # Uses Gauss-Chebyshev quadrature to integrate over all dimensions
  
  grid        = plan.grid
  multi_index = plan.multi_index
  domain      = plan.domain

  iim = smolyak_inverse_interpolation_matrix(grid,multi_index,domain) 

  d = size(grid,2)

  e = zeros(1,size(grid,1))
  e[1] = π^d
  w = e*iim

  y = zeros(size(grid,1))
  for i in eachindex(y)
    integrating_weights = sqrt(1.0-normalize_node(grid[i,1],domain[:,1])^2)
    for j = 2:d
      integrating_weights *= sqrt(1.0-normalize_node(grid[i,j],domain[:,j])^2)
    end
    y[i] = f(grid[i,:])*integrating_weights
  end

  scale_factor = (domain[1,1]-domain[2,1])/2
  for i in 2:d
    scale_factor = scale_factor*(domain[1,i]-domain[2,i])/2
  end

  return (w*y)[1] * scale_factor

end

########################################################################
########################################################################

"""
Uses the ```node_type``` function to construct the ```d```-dimensional full Smolyak grid with approximation layer
```mu``` and ```domain```.  If ```domain``` is not provided, then the approximation domain defaults to [1.0,-1.0]^d.
Returns the full Smolyak grid and the associated multi index.

Signature
=========

full_grid, mi = smolyak_grid_full(chebyshev_extrema,d,mu,domain)

Examples
========
```
julia> full_grid, mi = smolyak_grid_full(chebyshev_extrema,2,2)
julia> full_grid, mi = smolyak_grid_full(chebyshev_extrema,2,[2,1])
julia> full_grid, mi = smolyak_grid_full(chebyshev_extrema,2,2,[3.0 1.5; 2.0 0.5])
julia> full_grid, mi = smolyak_grid_full(chebyshev_extrema,2,[2,2],[3.0 1.5; 2.0 0.5])
```
"""
function smolyak_grid_full(node_type::Function,d::S,mu::S,domain=[ones(1,d);-ones(1,d)]) where {S<:Integer}

  T = typeof(1.0)

  multi_index        = generate_multi_index(d,mu)
  unique_multi_index = sort(unique(multi_index))
  unique_node_number = m_i.(unique_multi_index)

  # Create base nodes to be used in the sparse grid

  base_nodes   = Array{Array{T,1},1}(undef,length(unique_node_number))
  base_weights = Array{Array{T,1},1}(undef,length(unique_node_number))
  @inbounds for i in eachindex(unique_node_number)
    base_nodes[i] = node_type(unique_node_number[i])
  end

  # Select the relevant polynomials from the multi index

  ii = (sum(multi_index,dims=2) .>= max(d,mu+1)).*(sum(multi_index,dims=2) .<= d+mu)
  multi_index_full = zeros(S,sum(ii),size(multi_index,2))
  j = 1
  @inbounds for i in axes(multi_index,1)
    if ii[i] == true
      multi_index_full[j,:] = multi_index[i,:]
      j += 1
    end
  end

  # Construct the sparse grid from the nodes

  nodes = Array{T,2}(undef,determine_grid_size_full(multi_index_full))
  l = 1
  for j in axes(multi_index_full,1)
    new_nodes = base_nodes[multi_index_full[j,1]]  # Here new_nodes is a 1d array
    for i = 2:d
      new_nodes = combine_nodes(new_nodes,base_nodes[multi_index_full[j,i]])  # Here new_nodes becomes a 2d array
    end
    m = size(new_nodes,1)
    nodes[l:l+m-1,:] = new_nodes
    l += m
  end

  scale_nodes!(nodes,domain)

  return nodes, multi_index_full

end

function determine_grid_size_full(mi::Array{S,2}) where {S<:Integer} # Internal function, not exported

  temp = similar(mi)

  @inbounds for i in axes(mi,1)
    @inbounds for j in axes(mi,2)
      if mi[i,j] == 1
        temp[i,j] = 1
      else
        temp[i,j] = 2^(mi[i,j]-1)+1
      end
    end
  end

  s = 0
  @inbounds for i in axes(mi,1)
    t = 1
    @inbounds for j in axes(mi,2)
      t *= temp[i,j]
    end
    s += t
  end

  return (s, size(mi,2))

end

"""
Creates the master_index from a ```multi_index```.  Returns a matrix of integers.

Signature
=========

master_i = master_index(multi_index)

Example
=======
```
julia> mi = generate_multi_index(2,2)
julia> master_i = master_index(mi)
[ 1  1
  2  3
  5  5
 10  3
 13  9
 22  5]
```
"""
function master_index(multi_index::Array{S,2}) where {S<:Integer}

  temp_ind   = similar(multi_index)
  master_ind = zeros(S,size(multi_index,1),2)

  @inbounds for i in eachindex(multi_index)
    if multi_index[i] == 1
      temp_ind[i] = 1
    else
      temp_ind[i] = m_i(multi_index[i])
    end
  end

  master_ind[1,1] = 1
  master_ind[1,2] = prod(temp_ind[1,:])

  @inbounds for i = 2:size(master_ind,1)
    master_ind[i,1] = master_ind[i-1,1] + master_ind[i-1,2]
    master_ind[i,2] = prod(temp_ind[i,:])
  end

  return master_ind

end

function cheb_poly(order::S,x::R) where {S<:Integer,R<:Number} # Internal function, not exported

  p  = one(R)
  p1 = zero(R)
  p2 = zero(R)

  for i = 2:order+1
    if i == 2
      p1, p = p, x
    else
      p2, p1 = p1, p
        p  = 2*x*p1-p2
    end
  end

  return p

end

function prod_cjs(max_grid::Union{Array{T,1},Array{T,2}},min_grid::Union{Array{T,1},Array{T,2}},poly_grid::Array{T,2}) where {T<:AbstractFloat} # Internal function, not exported

  cjs = ones(size(poly_grid))

  @inbounds for i in axes(poly_grid,1)
    @inbounds for j in axes(poly_grid,2)
      if poly_grid[i,j] == max_grid[j] || poly_grid[i,j] == min_grid[j]
        cjs[i,j] *= 2
      end
    end
  end

  return prod(cjs,dims=2)

end

function compute_scale_factor(multi_index::Array{S,1}) where {S<:Integer} # Internal function, not exported

  scale_factor = 1.0

  @inbounds for i in eachindex(multi_index)
    if multi_index[i] > 1
      scale_factor *= 2.0/(m_i(multi_index[i]) - 1)
    end
  end

  return scale_factor

end

"""
Uses Chebyshev polynomials as basis functions to compute the weights in a full Smolyak polynomial approximation given the
approximation sample, ```y```, the approximation ```grid```, the ```multi_index```, and the approximation 
```domain``` (defaults to [1.0,-1.0]^d).  Returns a vector of vectors containing the weights in the Smolyak polynomial.

Signature
=========

w = smolyak_weights_full(y,grid,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid_full(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> f(x) = sum(x.^2)
julia> y = [f(g[i,:]) for i in axes(g,1)]
julia> w = smolyak_weights_full(y,g,mi,[1.0 1.0; 0.0 0.0])
[[-0.625, -0.5, -0.125]
 [0.625, 0.49999999999999994, 0.12499999999999992, -1.6653345369377348e-16, 2.7755575615628914e-17]
 [-0.625, -0.5, -0.125]
 [0.75, 0.5, 0.125, 0.5, 0.0, 0.0, 0.125, 0.0, 0.0]
 [0.625, 0.49999999999999994, 0.12499999999999992, -1.6653345369377348e-16, 2.7755575615628914e-17]]
```
"""
function smolyak_weights_full(y_f::Array{T,1},grid::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain=[ones(1,size(grid,2));-ones(1,size(grid,2))]) where {S<:Integer, T<:AbstractFloat}

  grid = copy(grid)
  for i in axes(grid,2)
    grid[:,i] = normalize_node(grid[:,i],domain[:,i])
  end

  mi = sum(multi_index,dims=2)
  d  = size(multi_index,2)
  mu = maximum(mi)-d

  max_grid = maximum(grid,dims=1)
  min_grid = minimum(grid,dims=1)

  weights = Array{Array{T,1},1}(undef,size(multi_index,1))
  g_ind = master_index(multi_index)

  @inbounds for i in axes(g_ind,1) # This loops over the number of polynomials
    ws = zeros(g_ind[i,2])
    poly_grid = grid[g_ind[i,1]:g_ind[i,1]+g_ind[i,2]-1,:]
    poly_y    = y_f[g_ind[i,1]:g_ind[i,1]+g_ind[i,2]-1]
      if size(grid,1) == 1 # This is to accommodate the mu = 0 case
        cjs_prod = ones(size(poly_grid,1))
      else
      cjs_prod = prod_cjs(max_grid,min_grid,poly_grid)
    end
    @inbounds for l = 1:g_ind[i,2] # This loops over the weights
      ll = CartesianIndices(Tuple(m_i.(multi_index[i,:])))[l]
      @inbounds for j = 1:g_ind[i,2] # This loops over the nodes
        rhs_term = cheb_poly(ll[1]-1,poly_grid[j,1])*poly_y[j]
        @inbounds for k = 2:size(poly_grid,2) # This loops over the polynomial terms in the product
          rhs_term *= cheb_poly(ll[k]-1,poly_grid[j,k])
        end
        ws[l] += rhs_term/cjs_prod[j]
      end
      scale_factor = compute_scale_factor(multi_index[i,:])
      ws[l] = scale_factor*(1/cjs_prod[l])*ws[l]
    end
    weights[i] = (-1)^(d+mu-mi[i])*factorial(d-1)/(factorial(d+mu-mi[i])*factorial(-1-mu+mi[i]))*ws
  end

  return weights

end

"""
Evaluates a full Smolyak polynomial formed using Chebyshev basis functions, given the ```weights```, the ```point``` at
which to evaluate the polynomial, the ```multi_index```, and the approximation ```domain``` (defaults to [1.0,-1.0]^d).
Returns a scalar.

Signature
=========

yhat = smolyak_evaluate_full(weights,point,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid_full(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = smolyak_weights_full(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> yhat = smolyak_evaluate_full(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0])
0.5953026581237829
```
"""
function smolyak_evaluate_full(weights::Array{Array{T,1},1},point::Array{R,1},multi_index::Array{S,2},domain=[ones(1,length(point));-ones(1,length(point))]) where {S<:Integer,R<:Number,T<:AbstractFloat}

  d = size(multi_index,2)
  point = copy(point)
  for i in eachindex(point)
    point[i] = normalize_node(point[i],domain[:,i])
  end

  evaluated_polynomials = zeros(size(multi_index,1))

  @inbounds for i in axes(multi_index,1) # This loops over the number of polynomials
    @inbounds for l in eachindex(weights[i])
      ll = CartesianIndices(Tuple(m_i.(multi_index[i:i,:])))[l]
      temp = weights[i][l]*cheb_poly(ll[1]-1,point[1])
      @inbounds for k = 2:d
        temp *= cheb_poly(ll[k]-1,point[k])
      end
      evaluated_polynomials[i] += temp
    end
    #evaluated_polynomials[i] *= (-1)^(d+mu-mi[i])*factorial(d-1)/(factorial(d+mu-mi[i])*factorial(-1-mu+mi[i]))
  end

  return sum(evaluated_polynomials)

end

function deriv_cheb_poly(order::S,x::R) where {S<:Integer,R<:Number} # Internal function, not exported

  p0 = one(R)
  p1 = zero(R)
  p2 = zero(R)
  pd0 = zero(R)
  pd1 = zero(R)
  pd2 = zero(R)

  for i = 2:order+1
    if i == 2
      p1, p0 = p0, x
      pd1, pd0 = pd0, one(R)
    else
      p2, p1 = p1, p0
      pd2, pd1 = pd1, pd0
        p0  = 2*x*p1-p2
      pd0 = 2*p1+2*x*pd1-pd2
    end
  end

  return pd0

end

function _smolyak_derivative_full(weights::Array{Array{T,1},1},point::Array{R,1},multi_index::Array{S,2},pos::S) where {S<:Integer,R<:Number,T<:AbstractFloat} # Internal function, not exported

  mi = sum(multi_index,dims=2)
  d  = size(multi_index,2)
  mu = maximum(mi)-d

  evaluated_polynomials = zeros(size(multi_index,1))

  for i in axes(multi_index,1) # This loops over the number of polynomials
    for l in eachindex(weights[i])
      ll = CartesianIndices(Tuple(m_i.(multi_index[i:i,:])))[l]
      if pos == 1
        temp = weights[i][l]*deriv_cheb_poly(ll[1]-1,point[1])
      else
        temp = weights[i][l]*cheb_poly(ll[1]-1,point[1])
      end
      for k = 2:d
        if k == pos
          temp *= deriv_cheb_poly(ll[k]-1,point[k])
        else
          temp *= cheb_poly(ll[k]-1,point[k])
        end
      end
      evaluated_polynomials[i] += temp
    end
    #evaluated_polynomials[i] *= (-1)^(d+mu-mi[i])*factorial(d-1)/(factorial(d+mu-mi[i])*factorial(-1-mu+mi[i]))
  end

  evaluated_derivative = sum(evaluated_polynomials)

  return evaluated_derivative

end

"""
Computes the partial derivative of a full Smolyak polynomial formed using Chebyshev basis functions, given the 
```weights```, the ```point``` at which to evaluate the polynomial, the ```multi_index```, the approximation 
```domain```, and the index of the variable to differentiate with respect to, ```pos```.
Returns a scalar.

Signature
=========

deriv = smolyak_derivative_full(weights,point,multi_index,domain,pos)

Example
=======
```
julia> g,mi = smolyak_grid_full(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = smolyak_weights_full(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> deriv1 = smolyak_derivative_full(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0],1)
0.4036816953853277
julia> deriv2 = smolyak_derivative_full(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0],2)
0.5250627466157655
```
"""
function smolyak_derivative_full(weights::Array{Array{T,1},1},point::Array{R,1},multi_index::Array{S,2},domain::Array{T,2},pos::S) where {S<:Integer,R<:Number,T<:AbstractFloat}

  d = size(multi_index,2)
  point = copy(point)
  for i = 1:d
    point[i] = normalize_node(point[i],domain[:,i])
  end

  evaluated_derivative = _smolyak_derivative_full(weights,point,multi_index,pos)*(2.0/(domain[1,pos]-domain[2,pos]))

  return evaluated_derivative

end

"""
Computes the gradient a full Smolyak polynomial formed using Chebyshev basis functions, given the ```weights```, the
```point``` at which to evaluate the polynomial, the ```multi_index```, and the approximation ```domain```.  
Returns a one-row matrix.

Signature
=========

grad = smolyak_gradient(weights,point,multi_index,domain)

Example
=======
```
julia> g,mi = smolyak_grid_full(chebyshev_extrema,2,2,[1.0 1.0; 0.0 0.0])
julia> y = g[:,1].^0.3.*g[:,2].^0.7
julia> w = smolyak_weights_full(y,g,mi,[1.0 1.0; 0.0 0.0])
julia> grad = smolyak_gradient_full(w,[0.37,0.71],mi,[1.0 1.0; 0.0 0.0])
[0.403682  0.525063]
```
"""
function smolyak_gradient_full(weights::Array{Array{T,1},1},point::Array{R,1},multi_index::Array{S,2},domain::Array{T,2}) where {S<:Integer,R<:Number,T<:AbstractFloat}

  d = size(multi_index,2)

  gradient = Array{R,2}(undef,1,d)
  for i = 1:d
    gradient[i] = smolyak_derivative_full(weights,point,multi_index,domain,i)
  end

  return gradient

end