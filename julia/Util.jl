module Util

using Random

export sample_linear_integer_weights

"""
Returns an index from 1:length(w) with probability proportional to w."

s must be precomputed to be sum(w).
"""
function sample_linear_integer_weights(rng::MersenneTwister, w::Vector{Int64}, s::Int64)
    i = rand(rng, 1:s)
    cs = 0
    for j = 1:(length(w) - 1)
        cs += w[j]
        if i <= cs
            return j
        end
    end
    length(w)
end

end
