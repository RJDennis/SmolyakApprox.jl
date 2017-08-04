function m_i{S<:Integer}(multi_index::Array{S,1})

  m_node_number = similar(multi_index)

  for i in eachindex(multi_index)

    if multi_index[i] == 1
      m_node_number[i] = 1
    else
      m_node_number[i] = 2^(multi_index[i]-1)+1
    end

  end

  return m_node_number

end
