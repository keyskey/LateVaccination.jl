# Run this file with
# julia -p n parallel_run.jl  n: number of process

@everywhere module Society
    using Statistics

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

    count_strategy_fraction(society::SocietyType) = length(filter(strategy -> strategy == "PV", society.strategy))/society.total_population 

    count_SAP(society::SocietyType) = Statistics.mean(society.point)
end

@everywhere module Epidemics
    using StatsBase
    using Random
    using Printf
    using ..Society

    function initialize_epidemics(society::Society.SocietyType, beta::Float64, gamma::Float64, delta::Float64, num_initial_i::Int)
        society = set_total_probability(reset_elapse_days(initialize_state(society, num_initial_i)), beta, gamma, delta)

        return society
    end
    
    function reset_elapse_days(society::Society.SocietyType)
        society.elapse_days = 0

        return society
    end

    function initialize_state(society::Society.SocietyType, num_initial_i::Int)
        society.survivors = Vector(1:society.total_population)  # Need to remove IM state

        # Pre-emptive Vは最初からIM stateなのでinitial iはLVから選ぶ
        LV_id = [id for (id, strategy) in enumerate(society.strategy) if strategy == "LV"]
        initial_i::Vector{Int} = StatsBase.self_avoid_sample!(LV_id, Vector(1:num_initial_i))
        society.num_s = society.num_im = society.num_i = society.num_r = 0  # Need to count!!
        for id in 1:society.total_population
            if id in initial_i
                society.state[id] = "I"
                society.num_i += 1
            elseif society.strategy[id] == "PV"
                society.state[id] = "IM"
                society.num_im += 1
                deleteat!(society.survivors, findfirst(isequal(id), society.survivors))
            else
                society.state[id] = "S"
                society.num_s += 1
            end
        end

        return society
    end
    
    function set_total_probability(society::Society.SocietyType, beta::Float64, gamma::Float64, delta::Float64)
        Ps_i = beta * society.num_i
        Pi_r = gamma
        P_vaccinate =  delta * society.num_i/(society.num_im + 1)
        society.total_probability = (Ps_i + P_vaccinate) * society.num_s + Pi_r * society.num_i

        return society
    end

    function one_season(society::Society.SocietyType, beta::Float64, gamma::Float64, delta::Float64, cr::Float64, season::Int)
        for timestep in 1:10000000
            rand_num::Float64 = rand()
            accum_probability::Float64 = 0.0
            someone_selected::Bool = false

            # State change loop
            for id in society.survivors
                if society.state[id] == "S"
                    accum_probability += beta * society.num_i/society.total_probability
                elseif society.state[id] == "I"
                    accum_probability += gamma/society.total_probability
                else
                    error("Error in state change, my state doesn't match with S/I")
                end

                if rand_num <= accum_probability
                    society = state_change(society, id)
                    someone_selected = true
                    break  # Escape from survivors loop
                end
            end

            # Vaccination loop, only health free riders are selected
            if someone_selected == false
                for id in society.survivors
                    if society.state[id] == "S" && society.strategy[id] == "LV"
                        accum_probability += delta * society.num_i/(society.num_im + 1)/society.total_probability
                        if rand_num <= accum_probability
                            society = vaccination(society, id)
                            someone_selected = true
                            break  # Escape from survivors loop
                        end
                    end
                end
            end

            # Error catch
            if someone_selected == false
                error("No one selected")
            end

            society = set_total_probability(society, beta, gamma, delta)
            #days = count_elapse_days(society)
            #fs, fim, fi, fr = Society.count_state_fraction(society)
            #@printf("Cr: %.1f Season: %i Step: %i Days: %.2f Fs: %.4f Fim: %.4f Fi: %.4f Fr: %.4f \n", cr, season, timestep, days, fs, fim, fi, fr)

            # Check conversion
            if society.num_i == 0
                break
            end
        end
        
        return society
    end

    function state_change(society::Society.SocietyType, id::Int)
        if society.state[id] == "S"
            society.state[id] = "I"
            society.num_s -= 1
            society.num_i += 1
        elseif society.state[id] == "I"
            society.state[id] = "R"
            society.num_i -= 1
            society.num_r += 1
            deleteat!(society.survivors, findfirst(isequal(id), society.survivors))
        else
            error("Error in desease spreading, R or the other state is picked up in gillespie method")
        end

        return society
    end

    function vaccination(society::Society.SocietyType, id::Int)
        society.state[id] = "IM"
        society.num_s -= 1
        society.num_im += 1
        deleteat!(society.survivors, findfirst(isequal(id), society.survivors))

        return society
    end

    function count_elapse_days(society::Society.SocietyType)
        society.elapse_days += log(1/rand())/society.total_probability

        return society.elapse_days
    end
end

@everywhere module Decision
    using StatsBase
    using ..Society

    function choose_initial_pv(total_population::Int, num_initial_pv::Int)
        initial_pv::Vector{Int} = StatsBase.self_avoid_sample!(Vector(1:total_population), Vector(1:num_initial_pv))

        return initial_pv
    end

    function initialize_strategy(society::Society.SocietyType, initial_pv::Vector{Int})
        for id in 1:society.total_population
            if id in initial_pv
                society.strategy[id] = "PV"
            else
                society.strategy[id] = "LV"
            end
        end

        return society
    end

    function count_payoff(society::Society.SocietyType, cr::Float64)
        for id in 1:society.total_population
            if society.state[id] == "IM"
                society.point[id] = -cr
            elseif society.strategy[id] == "LV"
                if society.state[id] == "S"
                    society.point[id] = 0
                elseif society.state[id] == "R"
                    society.point[id] = -1
                end
            end
        end

        return society
    end

    function update_strategy(society::Society.SocietyType)
        next_strategy::Vector{AbstractString} = copy(society.strategy)
        for id in 1:society.total_population
            opp_id = rand(1:society.total_population)
            while opp_id == id
                opp_id = rand(1:society.total_population)
            end
            
            if society.strategy[opp_id] != society.strategy[id] && rand() < 1/( 1 + exp( ( society.point[id] - society.point[opp_id] )/0.1 ) )
                next_strategy[id] = society.strategy[opp_id]
            end
        end
        
        society.strategy = copy(next_strategy)

        return society
    end
end

@everywhere module Simulation
    using Distributed
    using Random
    using CSV
    using DataFrames
    using Statistics
    using Printf
    using ..Epidemics
    using ..Decision
    using ..Society

    # Calculation for a pair of Cr and Effectiveness. Continue until fv reaches to equilibrium
    function season_loop(society::Society.SocietyType, beta::Float64, gamma::Float64, delta::Float64 , cr::Float64, num_initial_i::Int, initial_pv::Vector{Int})
        # Initialization
        society = Decision.initialize_strategy(society, initial_pv)
        fpv0 = Society.count_strategy_fraction(society)
        fpv_hist = [fpv0]

        # Season loop
        for season in 1:500
            # Initialization
            society = Epidemics.initialize_epidemics(society, beta, gamma, delta, num_initial_i)
            # fs0, fim0, fi0, fr0 = Society.count_state_fraction(society)
            fpv = Society.count_strategy_fraction(society)
            # @printf("Cr: %.1f Season: %i Step: 0 Days: %.2f Fs: %.4f Fim: %.4f Fpv: %.4f Fi: %.4f Fr: %.4f \n", cr, season, society.elapse_days, fs0, fim0, fpv, fi0, fr0)

            # One season
            society = Decision.update_strategy(Decision.count_payoff(Epidemics.one_season(society, beta, gamma, delta, cr, season), cr))
            global fs, Fim, fi, FES = Society.count_state_fraction(society)
            global SAP = Society.count_SAP(society)
            fpv_next = Society.count_strategy_fraction(society)
            push!(fpv_hist, fpv_next)

            # Check conversion
            if fpv_next * society.total_population >= society.total_population - num_initial_i  # Can't choose initial infected agent!
                global Fpv = fpv
                break
            elseif season >= 100 && Statistics.mean(fpv_hist[season-99:season]) - fpv_next < 0.001
                global Fpv = Statistics.mean(fpv_hist[season-98:season+1])
                break
            end
        end

        @printf("δ: %.1f Cr: %.1f Finished with FES: %.4f Fim(VC): %.4f Fpv: %.4f SAP: %.3f \n", delta, cr, FES, Fim, Fpv, SAP)

        return Dict("Delta" => delta, "Cr" => cr, "FES" => FES, "Fim" => Fim, "Fpv" => Fpv, "SAP" => SAP)
    end

    # Get data for one Cr-Effectiveness phase diagram
    function one_episode(total_population::Int, episode::Int)
        Random.seed!()
        beta::Float64 = 0.000086 
        gamma::Float64 = 1/3
        num_initial_i::Int = 5
        num_initial_pv::Int = div(total_population, 10)
        initial_pv::Vector{Int} = Decision.choose_initial_pv(total_population, num_initial_pv)
        society = Society.SocietyType(total_population)

        for delta::Float64 in [0.2, 0.5, 0.8, 1.0]
            DataFrame(Cr = [], FES = [], Fim = [], Fpv = [], SAP = []) |> CSV.write("result$(episode)_delta_$(delta).csv")
            results::Vector{Dict{String, Float64}} = Distributed.pmap(cr -> Simulation.season_loop(society, beta, gamma, delta, cr, num_initial_i, initial_pv), 0:0.1:1.0)
            for result in results
                DataFrame(Cr = [result["Cr"]], FES = [result["FES"]], Fim = [result["Fim"]], Fpv = [result["Fpv"]], SAP = [result["SAP"]]) |> CSV.write("result$(episode)_delta_$(result["Delta"]).csv", append=true)
            end
        end
    end
end

using .Simulation

const num_episode = 1
const total_population = 10000

for episode = 1:num_episode
    Simulation.one_episode(total_population, episode)
end
