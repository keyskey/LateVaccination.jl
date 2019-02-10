module Society
    using Statistics
    export SocietyType

    mutable struct  SocietyType
        total_population::Int
        num_s::Int   # Free rider
        num_im::Int  # Immuned state
        num_i::Int
        num_r::Int
        state::Vector{AbstractString}      # S, IM, I, R
        strategy::Vector{AbstractString}   # PV or LV
        point::Vector{Float64}
        survivors::Vector{Int}
        total_probability::Float64
        elapse_days::Float64
        
        # Constructor
        SocietyType(total_population::Int) = new(
            total_population,                   # total_population
            0,                                  # num_s
            0,                                  # num_im
            0,                                  # num_i 
            0,                                  # num_r
            fill(" ", total_population),        # state
            fill(" ", total_population),        # strategy
            zeros(total_population),            # point
            Vector(1:total_population),         # survivors
            0,                                  # total_probability
            0                                   # elapse_days
            )
    end

    function count_state_fraction(society::SocietyType)
        n_tot = society.total_population
        fs = society.num_s/n_tot
        fim = society.num_im/n_tot
        fi = society.num_i/n_tot
        fr = society.num_r/n_tot
        
        return fs, fim, fi, fr
    end

    function count_strategy_fraction(society::SocietyType)
        fpv = length(filter(strategy -> strategy == "PV", society.strategy))/society.total_population
        flv = length(filter(strategy -> strategy == "LV", society.strategy))/society.total_population

        return fpv, flv
    end

    count_SAP(society::SocietyType) = Statistics.mean(society.point)
end
