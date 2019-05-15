function m_i(multi_index::Union{S,Array{S,1},Array{S,2}}) where {S<:Integer}

  if size(multi_index,2) != 1
	error("the multi index should be a vector")
  end

  m_node_number = copy(multi_index)

  if typeof(multi_index) == S
    if multi_index == 1
	  m_node_number = 1
    else
	  m_node_number = 2^(multi_index-1)+1
    end
  else
    for i = 1:length(multi_index)
      if multi_index[i] == 1
        m_node_number[i] = 1
      else
        m_node_number[i] = 2^(multi_index[i]-1)+1
      end
    end
  end

  return m_node_number

end
