# Run this file with
# julia -p n parallel_run.jl  n: number of process

#@everywhere module Society
module Society
    using Statistics
    export SocietyType, Params

    mutable struct  SocietyType
        num_tot::Int
        num_s::Int
        num_im::Int
        num_i::Int
        num_r::Int
        state::Vector{AbstractString}       # S, IM, I, R
        strategy::Vector{AbstractString}    # PV or LV
        point::Vector{Float64}
        survivors::Vector{Int}
        total_probability::Float64
        elapse_days::Float64

        # Constructor
        SocietyType(population::Int) = new(
            population,                         # num_tot
            0,                                  # num_s
            0,                                  # num_im
            0,                                  # num_i
            0,                                  # num_r
            fill(" ", population),              # state
            fill(" ", population),              # strategy
            zeros(population),                  # point
            Vector(1:population),               # survivors
            0,                                  # total_probability
            0,                                  # elapse_days
            )
    end

    struct Params
        beta::Float64
        gamma::Float64
        Cr::Float64
        effectiveness::Float64
        num_initial_i::Int
        initial_pv::Vector{Int}
    end

    function count_state_fraction(society::SocietyType)
        fs = society.num_s/society.num_tot
        fim = society.num_im/society.num_tot
        fi = society.num_i/society.num_tot
        fr = society.num_r/society.num_tot

        return fs, fim, fi, fr
    end

    count_fv(society::SocietyType) = length(filter(strategy -> strategy == "V", society.strategy))/society.num_tot

    count_SAP(society::SocietyType) = Statistics.mean(society.point)
end

#@everywhere module Epidemics
module Epidemics
    using StatsBase
    using Random
    using Printf
    using ..Society

    function initialize_epidemics(society::Society.SocietyType, params::Society.Params)
        society = initialize_state(society, params) |> (initialized_society -> set_total_probability(initialized_society, params))

        return society
    end
    
    function initialize_state(society::Society.SocietyType, params::Society.Params)
        society.survivors = Vector(1:society.num_tot)  # Need to remove vaccinated agents
        society.elapse_days = 0

        # Choose initial infected agent from NV
        NV_id = [id for (id, strategy) in enumerate(society.strategy) if strategy == "NV"]
        initial_i::Vector{Int} = StatsBase.self_avoid_sample!(NV_id, Vector(1:params.num_initial_i))
        society.num_s  = 0
        society.num_im = 0
        society.num_i  = 0
        society.num_r  = 0  
        
        for id in 1:society.num_tot
            if id in initial_i
                society.state[id] = "I"
                society.num_i += 1
            elseif society.strategy[id] == "V"
                if rand() < params.effectiveness
                    society.state[id] = "IM"
                    society.num_im += 1
                    deleteat!(society.survivors, findfirst(isequal(id), society.survivors))
                else
                    society.state[id] = "S"
                    society.num_s += 1
                end
            else
                society.state[id] = "S"
                society.num_s += 1
            end
        end

        return society
    end
    
    function set_total_probability(society::Society.SocietyType, params::Society.Params)
        Ps_i = params.beta * society.num_i
        Pi_r = params.gamma

        society.total_probability = Ps_i * society.num_s + Pi_r * society.num_i

        return society
    end

    function one_season(society::Society.SocietyType, params::Society.Params, season::Int)
        for timestep in 1:1000000
            rand_num = rand()
            accum_probability = 0.0
            someone_selected = false
            Ps_i = params.beta * society.num_i
            Pi_r = params.gamma

            for id in society.survivors
                if society.state[id] == "S"
                    accum_probability += Ps_i
                elseif society.state[id] == "I"
                    accum_probability += Pi_r
                end

                if rand_num <= accum_probability/society.total_probability
                    society = state_change(society, id)
                    someone_selected = true
                    break
                end
            end

            society = set_total_probability(society, params)
            days = count_elapse_days(society)
            fs, fim, fi, fr = Society.count_state_fraction(society)
            @printf("Cr: %.1f Effectiveness: %.1f Season: %i Step: %i Days: %.2f Fs: %.4f Fim: %.4f Fi: %.4f Fr: %.4f \n", 
                    params.Cr,
                    params.effectiveness,
                    season,
                    timestep,
                    days,
                    fs,
                    fim,
                    fi,
                    fr)
                
            # Error catch
            if someone_selected == false
                error("Error in one_season, no one selected !!")
            end
                
            # Check conversion
            if society.num_i == 0
                break
            end
        end
        
        return society
    end

    function state_change(society::Society.SocietyType, id::Int)
        if society.state[id] == "S"
            next_state = "I"
            society.num_s -= 1
            society.num_i += 1
        elseif society.state[id] == "I"
            next_state = "R"
            society.num_i -= 1
            society.num_r += 1
            deleteat!(society.survivors, findfirst(isequal(id), society.survivors))
        end
        society.state[id] = next_state

        return society
    end

    function count_elapse_days(society::Society.SocietyType)
        society.elapse_days += log(1/rand())/society.total_probability

        return society.elapse_days
    end
end

#@everywhere module Decision
module Decision
    using StatsBase
    using ..Society

    function choose_initial_v(population::Int, num_initial_v::Int)
        initial_v::Vector{Int} = StatsBase.self_avoid_sample!(Vector(1:population), Vector(1:num_initial_v))

        return initial_v
    end

    function initialize_strategy(society::Society.SocietyType, initial_v::Vector{Int})
        for id in 1:society.num_tot
            if id in initial_v
                society.strategy[id] = "V"
            else
                society.strategy[id] = "NV"
            end
        end

        return society
    end

    function count_payoff(society::Society.SocietyType, Cr::Float64)
        for id in 1:society.num_tot
            if society.strategy[id] == "V"
                # Healthy vaccinators
                if society.state[id] in ("S", "IM") 
                    society.point[id] = -Cr
                # Infected vaccinators
                elseif society.state[id] == "R"
                    society.point[id] = -Cr-1
                end
            elseif society.strategy[id] == "NV"
                # Healthy non-vaccinators
                if society.state[id] == "S"
                    society.point[id] = 0
                # Infected non-vaccinators
                elseif society.state[id] == "R"
                    society.point[id] = -1
                end
            end
        end

        return society
    end

    function update_strategy(society::Society.SocietyType)
        next_strategy::Vector{AbstractString} = copy(society.strategy)
        for id in 1:society.num_tot
            opp_id = rand(1:society.num_tot)
            while opp_id == id
                opp_id = rand(1:society.num_tot)
            end
            
            if society.strategy[opp_id] != society.strategy[id] && rand() < 1/( 1 + exp( ( society.point[id] - society.point[opp_id] )/0.1 ) )
                next_strategy[id] = society.strategy[opp_id]
            end
        end
        
        society.strategy = copy(next_strategy)

        return society
    end
end

#@everywhere module Simulation
module Simulation
    using Distributed
    using Random
    using CSV
    using DataFrames
    using Statistics
    using Printf
    using ..Epidemics
    using ..Decision
    using ..Society

    # Calculation for a pair of (Cr, Effectiveness).
    function season_loop(society::Society.SocietyType, params::Society.Params)
        ###################### Initialization #####################
        society = Decision.initialize_strategy(society, params.initial_pv)
        fv0 = Society.count_fv(society)
        fes_hist = []
        sap_hist = []
        fv_hist = [fv0]
        season = 0
        ###########################################################

        while true 
            ######################## Initialization ########################
            season += 1
            society = Epidemics.initialize_epidemics(society, params)
            fs0, fim0, fi0, fr0 = Society.count_state_fraction(society)
            fv_this_season = Society.count_fv(society)
            @printf("Cr: %.1f Effectiveness: %.1f Season: %i Step: %i Days: %.1f Fs: %.4f Fim: %.4f Fi: %.4f Fr: %.4f \n", 
                    params.Cr,
                    params.effectiveness,
                    season,
                    0,
                    society.elapse_days,
                    fs0,
                    fim0,
                    fi0,
                    fr0)
            ################################################################

            # One season
            society = Epidemics.one_season(society, params, season) |> (result -> Decision.count_payoff(result, params.Cr)) |> Decision.update_strategy
            fs, fim, fi, fes = Society.count_state_fraction(society)
            sap = Society.count_SAP(society)
            fv_next_season = Society.count_fv(society)
            push!(fes_hist, fes)
            push!(sap_hist, sap)
            push!(fv_hist, fv_next_season)
            
            # Check conversion
            if (1 - fv_next_season) * society.num_tot <= params.num_initial_i  # If number of LVs is less than num_initial_i, you can't choose initial infected agent!
                # Use the value at the previous season as a final solution
                global FES = fes_hist[end]
                global SAP = sap_hist[end]
                global Fv  = fv_hist[end-1]
                break

            elseif (season >= 100 && abs(Statistics.mean(fv_hist[season-99:season]) - fv_next_season) < 0.001) || season == 500
                global FES = Statistics.mean(fes_hist[season-99:season])
                global SAP = Statistics.mean(sap_hist[season-99:season])
                global Fv  = Statistics.mean(fv_hist[season-99:season])
                break
            end
        end

        @printf("Cr: %.1f Effectiveness: %.1f Finished with FES: %.4f Fv: %.4f SAP: %.3f \n", 
                params.Cr, 
                params.effectiveness,
                FES,
                Fv,
                SAP)

        result = Dict("Cr"            => params.Cr,
                      "effectiveness" => params.effectiveness,
                      "FES"           => FES,
                      "Fv"            => Fv,
                      "SAP"           => SAP)

        return result
    end

    # Get data for one Cr-Effectiveness phase diagram
    function one_episode(population::Int, episode::Int)
        Random.seed!()
        beta::Float64 = 0.000086 
        gamma::Float64 = 1/3
        num_initial_i::Int = 5
        num_initial_v::Int = div(population, 2)
        initial_v::Vector{Int} = Decision.choose_initial_v(population, num_initial_v)
        society = Society.SocietyType(population)

        DataFrame(Cr = [], effectiveness = [], FES = [], Fv = [], SAP = []) |> CSV.write("conventional_result$(episode).csv")
        
        for effectiveness in 0:0.1:1.0
            results::Vector{Dict{String, Float64}} = Distributed.pmap(Cr -> Simulation.season_loop(society, Params(beta, gamma, Cr, effectiveness, num_initial_i, initial_v)), 0:0.1:1.0)
            #results::Vector{Dict{String, Float64}} = map(Cr -> Simulation.season_loop(society, Params(beta, gamma, Cr, effectiveness, num_initial_i, initial_v)), 0:0.1:0.1)
    
            for result in results
                DataFrame(Cr            = [result["Cr"]], 
                          effectiveness = [result["effectiveness"]], 
                          FES           = [result["FES"]],
                          Fv            = [result["Fv"]], 
                          SAP           = [result["SAP"]]) |> CSV.write("conventional_result$(episode).csv", append=true)
            end
        end
    end
end

using .Simulation

const num_episode = 1
const population = 10000

for episode = 1:num_episode
    Simulation.one_episode(population, episode)
end
