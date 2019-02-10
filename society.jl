module Society
    using Statistics
    using Random
    export SocietyType

    mutable struct  SocietyType
        total_population::Int
        num_s::Int   # Failed vaccinators + Non vaccinators
        num_im::Int  # Immuned state (= Successful vaccinators)
        num_i::Int
        num_r::Int
        num_v::Int
        num_nv::Int
        num_sv::Int
        state::Vector{AbstractString}      # S, IM, E, I, R
        strategy::Vector{AbstractString}   # V or NV
        point::Vector{Float64}
        survivors::Vector{Int}
        opponents::Vector{Int}
        total_probability::Float64
        elapse_days::Float64
        
        # Constructor
        SocietyType(total_population::Int) = new(
            total_population,                   # total_population
            0,                                  # num_s
            0,                                  # num_im
            0,                                  # num_i 
            0,                                  # num_r
            0,                                  # num_v
            0,                                  # num_nv
            0,                                  # num_sv
            fill(" ", total_population),        # state
            fill(" ", total_population),        # strategy
            zeros(total_population),            # point
            Vector(1:total_population),         # survivors
            Random.shuffle(1:total_population), # opponents
            0,                                  # total_probability
            0                                   # elapse_days
            )
    end

    function count_fraction(society)
        n_tot = society.total_population
        fs = society.num_s/n_tot
        fim = society.num_im/n_tot
        fi = society.num_i/n_tot
        fr = society.num_r/n_tot
        fv =  society.num_v/n_tot
        fnv = society.num_nv/n_tot
        
        fs, fim, fi, fr, fv, fnv
    end

    count_SAP(society) = Statistics.mean(society.point)
end