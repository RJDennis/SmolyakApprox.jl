function m_i(multi_index::S) where {S<:Integer}

  if multi_index == 1
    m_node_number = one(S)
  else
    m_node_number = 2^(multi_index-1)+1
   end

  return m_node_number

end

function m_i(multi_index::Array{S,1}) where {S<:Integer}

  m_node_number = Array{S,1}(undef,length(multi_index))

  @inbounds for i in eachindex(multi_index)
    m_node_number[i] = m_i(multi_index[i])
  end

  return m_node_number

end

function m_i(multi_index::Array{S,2}) where {S<:Integer}

  m_node_number = Array{S,2}(undef,size(multi_index))

  @inbounds for i in eachindex(multi_index)
    m_node_number[i] = m_i(multi_index[i])
  end

  return m_node_number

end