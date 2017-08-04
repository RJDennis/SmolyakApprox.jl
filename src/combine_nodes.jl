function combine_nodes(nodes1,nodes2)  # nodes1 can be a 1d or 2d array; nodes2 is a 1d array

  combined_nodes = [kron(ones(size(nodes2,1)),nodes1) kron(nodes2,ones(size(nodes1,1)))]

  return combined_nodes

end
