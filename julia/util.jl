using Random
using Distributions
using StatsBase
using DelimitedFiles
using Dates
using JSON2


############## structures.jl should be loaded first ############

JSON2.@format Parameters noargs


function load_parameters_from_json(filename) :: Parameters
    str = open(f -> read(f, String), filename)
    params = JSON2.read(str, Parameters)
    validate(params)
    params
end


### VALIDATION TO ENSURE CODE CORRECTNESS IN THE FACE OF OPTIMIZATIONS ###

function validate(p::Parameters)
    @assert p.t_final !== nothing
    @assert p.t_output !== nothing

    @assert p.n_bstrains !== nothing
    @assert p.n_hosts_per_bstrain !== nothing
    @assert p.n_vstrains !== nothing
    @assert p.n_particles_per_vstrain !== nothing
    @assert p.n_protospacers !== nothing

    @assert p.u_n_spacers_max !== nothing
    @assert p.p_crispr_failure_prob !== nothing
    @assert p.q_spacer_acquisition_prob !== nothing
    @assert p.r_growth_rate !== nothing
    @assert p.K_carrying_capacity !== nothing
    @assert p.beta_burst_size !== nothing
    @assert p.phi_adsorption_rate !== nothing
    @assert p.m_viral_decay_rate !== nothing
    @assert p.mu_viral_mutation_rate !== nothing
    @assert p.rho_c_density_cutoff !== nothing
    @assert p.d_death_rate !== nothing
    @assert p.g_immigration_rate !== nothing
end

function validate(s::State)
    validate(s.bstrains)
    validate(s.vstrains)
end

function validate(strains::Strains)
    @assert strains.total_abundance == sum(strains.abundance)
    @assert strains.next_id > maximum(strains.ids)
    @assert length(strains.abundance) == length(strains.ids)
    @assert length(strains.spacers) == length(strains.ids)
end


### RNG UTILITY FUNCTIONS ##

"""
Returns an index from 1:length(w) with probability proportional to w.

s must be precomputed to be sum(w).
"""
function sample_linear_integer_weights(rng::MersenneTwister, w::Vector{UInt64}, s::UInt64)
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

"""
Removes an item in the middle of an array that does not need to be kept ordered in constant time.

The item is replaced with the item at the end of the array, and then the item at the end of the
array is removed.
"""
function swap_with_end_and_remove!(a, index)
    if index != lastindex(a)
        setindex!(a, a[lastindex(a)], index)
    end
    pop!(a)
    nothing
end

function remove_strain!(strains, index)
    # This is only used when a strain has gone extinct
    @assert strains.abundance[index] == 0

    @debug "Removing strain" id=strains.ids[index] index=index

    swap_with_end_and_remove!(strains.ids, index)
    swap_with_end_and_remove!(strains.abundance, index)
    swap_with_end_and_remove!(strains.spacers, index)
end
