# This function relates to the ansiotropic case

function generate_multi_index(d::S,mu::Array{S,1}) where {S<:Integer}

  max_mu = maximum(mu)
  multi_index = ones(S,1,d)
  for i = 2:(max_mu+1)^d

    candidate_index = collect(ind2sub((repeat([max_mu+1],inner = d)...),i))
    if sum(candidate_index) <= d+max_mu && sum(candidate_index .<= mu+1) == d
      multi_index = [multi_index; candidate_index']
    end

  end

  return multi_index

end

# The function below relates to the isotropic case

function generate_multi_index(d::S,mu::S) where {S<:Integer}

  multi_index = ones(S,1,d)
  for i = 2:(mu+1)^d

    candidate_index = collect(ind2sub((repeat([mu+1],inner = d)...),i))
    if sum(candidate_index) <= d+mu
      multi_index = [multi_index; candidate_index']
    end

  end

  return multi_index

end
