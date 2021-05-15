function combine_nodes(nodes1::Union{Array{T,1},Array{T,2}},nodes2::Array{T,1}) where {T<:AbstractFloat}  # nodes1 can be a 1d or 2d array; nodes2 is a 1d array

  n1 = size(nodes1,1)
  n2 = size(nodes1,2)
  n3 = length(nodes2)

  combined_nodes = Array{T,2}(undef,n1*n3,n2+1)

  for i = 1:n3
    combined_nodes[(i-1)*n1+1:i*n1,1:n2] = nodes1
  end
  for i = 1:n1
    for j = 1:n3
      combined_nodes[(j-1)*n1+i,n2+1] = nodes2[j]
    end
  end

  return combined_nodes

end

#function combine_nodes(nodes1::Union{Array{T,1},Array{T,2}},nodes2::Array{T,1}) where {T<:AbstractFloat}  # nodes1 can be a 1d or 2d array; nodes2 is a 1d array

#  combined_nodes = [kron(ones(size(nodes2,1)),nodes1) kron(nodes2,ones(size(nodes1,1)))]

#  return combined_nodes

#end
