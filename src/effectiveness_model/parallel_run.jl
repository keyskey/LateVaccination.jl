# Run this file with
# julia -p n parallel_run.jl  n: number of process

#@everywhere module Society
module Society
    using Statistics

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
        vaccinated::Vector{Bool}

        # Constructor
        SocietyType(population::Int) = new(
            population,                         # num_tot
            0,                                  # num_s
            0,                                  # num_im
            0,                                  # num_i
            0,                                  # num_r
            fill(" ", population),        # state
            fill(" ", population),        # strategy
            zeros(population),            # point
            Vector(1:population),         # survivors
            0,                                  # total_probability
            0,                                  # elapse_days
            fill(false, population)       # vaccinated
            )
    end

    function count_state_fraction(society::SocietyType)
        fs = society.num_s/society.num_tot
        fim = society.num_im/society.num_tot
        fi = society.num_i/society.num_tot
        fr = society.num_r/society.num_tot

        return fs, fim, fi, fr
    end

    count_fpv(society::SocietyType) = length(filter(strategy -> strategy == "PV", society.strategy))/society.num_tot

    count_SAP(society::SocietyType) = Statistics.mean(society.point)
end

#@everywhere module Epidemics
module Epidemics
    using StatsBase
    using Random
    using Printf
    using ..Society

    function initialize_epidemics(society::Society.SocietyType, beta::Float64, gamma::Float64, delta::Float64, effectiveness::Float64, num_initial_i::Int)
        initialized_society = initialize_state(society, num_initial_i, effectiveness) |> reset_elapse_days |> (result -> set_total_probability(result, beta, gamma, delta))

        return initialized_society
    end
    
    function reset_elapse_days(society::Society.SocietyType)
        society.elapse_days = 0

        return society
    end

    function initialize_state(society::Society.SocietyType, num_initial_i::Int, effectiveness::Float64)
        society.survivors = Vector(1:society.num_tot)  # Need to remove vaccinated agents
        society.vaccinated = fill(false, society.num_tot)

        # Choose initial infected agent from LV
        LV_id = [id for (id, strategy) in enumerate(society.strategy) if strategy == "LV"]
        initial_i::Vector{Int} = StatsBase.self_avoid_sample!(LV_id, Vector(1:num_initial_i))
        society.num_s = society.num_im = society.num_i = society.num_r = 0  # Need to count!!
        for id in 1:society.num_tot
            if id in initial_i
                society.state[id] = "I"
                society.num_i += 1
            elseif society.strategy[id] == "PV"
                society.vaccinated[id] = true
                if rand() < effectiveness
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
    
    function set_total_probability(society::Society.SocietyType, beta::Float64, gamma::Float64, delta::Float64)
        Ps_i = beta * society.num_i
        Pi_r = gamma
        P_vaccinate =  count_vaccination_probability(society, delta)
        num_healthy_non_v = count_num_healthy_non_v(society)

        society.total_probability = Ps_i * society.num_s + Pi_r * society.num_i + P_vaccinate * num_healthy_non_v

        return society
    end

    # Health(susceptible) agent should be included in survivors
    count_num_healthy_non_v(society::Society.SocietyType) = length([id for id in society.survivors if society.state[id] == "S" && society.vaccinated[id] == false])

    # P_vaccinate = δI/(V+1)
    count_vaccination_probability(society::Society.SocietyType, delta::Float64) = delta * society.num_i/(sum(society.vaccinated) + 1)

    function one_season(society::Society.SocietyType, beta::Float64, gamma::Float64, delta::Float64, cr::Float64, effectiveness::Float64, season::Int)
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
                end

                if rand_num <= accum_probability
                    society = state_change(society, id)
                    someone_selected = true
                    break  # Escape from survivors loop
                end
            end

            # Vaccination loop, only healthy non-vaccinators are selected
            if someone_selected == false
                for id in society.survivors
                    if society.state[id] == "S" && society.vaccinated[id] == false
                        P_vaccinate = count_vaccination_probability(society, delta)
                        accum_probability += P_vaccinate/society.total_probability
                        if rand_num <= accum_probability
                            society = vaccination(society, id, effectiveness)
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
            days = count_elapse_days(society)
            fs, fim, fi, fr = Society.count_state_fraction(society)
            @printf("Cr: %.1f Effectiveness: %.1f Season: %i Step: %i Days: %.2f Fs: %.4f Fim: %.4f Fi: %.4f Fr: %.4f \n", cr, effectiveness, season, timestep, days, fs, fim, fi, fr)

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
        end

        return society
    end

    function vaccination(society::Society.SocietyType, id::Int, effectiveness::Float64)
        society.vaccinated[id] = true
        if rand() < effectiveness
            society.state[id] = "IM"
            society.num_s -= 1
            society.num_im += 1
            deleteat!(society.survivors, findfirst(isequal(id), society.survivors))
        #= else
            state = "S"
            no change in num_s/im
        =#
        end

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

    function choose_initial_pv(population::Int, num_initial_pv::Int)
        initial_pv::Vector{Int} = StatsBase.self_avoid_sample!(Vector(1:population), Vector(1:num_initial_pv))

        return initial_pv
    end

    function initialize_strategy(society::Society.SocietyType, initial_pv::Vector{Int})
        for id in 1:society.num_tot
            if id in initial_pv
                society.strategy[id] = "PV"
            else
                society.strategy[id] = "LV"
            end
        end

        return society
    end

    function count_payoff(society::Society.SocietyType, Cr::Float64)
        for id in 1:society.num_tot
            if society.strategy[id] == "PV"
                # Healthy Pre-emptive Vaccinators
                if society.state[id] in ("S", "IM") 
                    society.point[id] = -Cr
                # Infected Pre-emptive Vaccinators
                elseif society.state[id] == "R"
                    society.point[id] = -Cr-1
                end
            elseif society.strategy[id] == "LV"
                # Healthy Late Vaccinators
                if society.state[id] == "IM" || (society.state[id] == "S" && society.vaccinated[id] == true)
                    society.point[id] = -Cr
                # Healthy Non-Vaccinators
                elseif society.state[id] == "S" && society.vaccinated[id] == false
                    society.point[id] = 0
                # Infected Late Vaccinators
                elseif society.state[id] == "R" && society.vaccinated[id] == true
                    society.point[id] = -Cr-1
                # Infected non-vaccinators
                elseif society.state[id] == "R" && society.vaccinated[id] == false
                    society.point[id] = -1
                end
            end
        end

        return society
    end

    function count_payoff(society::Society.SocietyType, Cpv::Float64, Cadd_rate::Float64)
        Clv::Float64 = Cpv*Cadd_rate

        for id in 1:society.num_tot
            if society.strategy[id] == "PV"
                # Healthy Pre-emptive Vaccinators
                if society.state[id] in ("S", "IM") 
                    society.point[id] = -Cpv
                # Infected Pre-emptive Vaccinators
                elseif society.state[id] == "R"
                    society.point[id] = -Cpv-1
                end
            elseif society.strategy[id] == "LV"
                # Healthy Late Vaccinators
                if society.state[id] == "IM" || (society.state[id] == "S" && society.vaccinated[id] == true)
                    society.point[id] = -Clv
                # Healthy Non-Vaccinators
                elseif society.state[id] == "S" && society.vaccinated[id] == false
                    society.point[id] = 0
                # Infected Late Vaccinators
                elseif society.state[id] == "R" && society.vaccinated[id] == true
                    society.point[id] = -Clv-1
                # Infected non-vaccinators
                elseif society.state[id] == "R" && society.vaccinated[id] == false
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

    # Calculation for a pair of (Cpv, Cadd_rate, Effectiveness).
    function season_loop(society::Society.SocietyType, beta::Float64, gamma::Float64, delta::Float64, Cpv::Float64, Cadd_rate::Float64, effectiveness::Float64, num_initial_i::Int, initial_pv::Vector{Int})
        # Initialization
        society = Decision.initialize_strategy(society, initial_pv)
        fpv0 = Society.count_fpv(society)
        fpv_hist = [fpv0]

        # Season loop
        season = 0
        while true
            # Initialization
            season += 1
            society = Epidemics.initialize_epidemics(society, beta, gamma, delta, effectiveness, num_initial_i)
            fs0, fim0, fi0, fr0 = Society.count_state_fraction(society)
            fpv = Society.count_fpv(society)
            @printf("Cpv: %.1f Cadd_rate: %.1f Season: %i Step: 0 Days: %.2f Fs: %.4f Fim: %.4f Fpv: %.4f Fi: %.4f Fr: %.4f \n", Cpv, Cadd_rate, season, society.elapse_days, fs0, fim0, fpv, fi0, fr0)

            # One season
            society = Epidemics.one_season(society, beta, gamma, delta, Cpv, effectiveness, season) |> (result -> Decision.count_payoff(result, Cpv, Cadd_rate)) |> Decision.update_strategy
            global fs, Fim, fi, FES = Society.count_state_fraction(society)
            global SAP = Society.count_SAP(society)
            fpv_next = Society.count_fpv(society)
            push!(fpv_hist, fpv_next)

            # Check conversion
            if fpv_next * society.num_tot >= society.num_tot - num_initial_i  # Can't choose initial infected agent!
                global Fpv = fpv
                break
            elseif season >= 100 && Statistics.mean(fpv_hist[season-99:season]) - fpv_next < 0.001
                global Fpv = Statistics.mean(fpv_hist[season-98:season+1])
                break
            end
        end

        @printf("δ: %.1f Cr: %.1f Cadd_rate: %.1f Effectiveness: %.1f Finished with FES: %.4f Fim(VC): %.4f Fpv: %.4f SAP: %.3f \n", delta, Cpv, Cadd_rate, effectiveness, FES, Fim, Fpv, SAP)

        return Dict("Delta" => delta, "Cr" => Cpv, "Cadd_rate" => Cadd_rate, "Effectiveness" => effectiveness, "FES" => FES, "Fim" => Fim, "Fpv" => Fpv, "SAP" => SAP)
    end

    # Get data for one Cr-Effectiveness phase diagram
    function one_episode(population::Int, episode::Int)
        Random.seed!()
        beta::Float64 = 0.000086 
        gamma::Float64 = 1/3
        num_initial_i::Int = 5
        num_initial_pv::Int = div(population, 10)
        initial_pv::Vector{Int} = Decision.choose_initial_pv(population, num_initial_pv)
        society = Society.SocietyType(population)

        #DataFrame(Delta = [], Cr = [], effectiveness = [], FES = [], Fim = [], Fpv = [], SAP = []) |> CSV.write("imperfect_result$(episode).csv")
        DataFrame(Delta = [], Cpv = [], Cadd_rate = [], effectiveness = [], FES = [], Fim = [], Fpv = [], SAP = []) |> CSV.write("imperfect_result$(episode).csv")

        for Cadd_rate::Float64 in 1.0:0.1:2.0
            for delta::Float64 in [0.2, 0.5, 1.0]
                for effectiveness in [1.0]  #0:0.1:1.0
                    #results::Vector{Dict{String, Float64}} = Distributed.pmap(Cpv -> Simulation.season_loop(society, beta, gamma, delta, Cpv, Cadd_rate, effectiveness, num_initial_i, initial_pv), 0:0.1:1.0)
                    results::Vector{Dict{String, Float64}} = map(Cpv -> Simulation.season_loop(society, beta, gamma, delta, Cpv, Cadd_rate, effectiveness, num_initial_i, initial_pv), 0:0.1:1.0)
                    for result in results
                        DataFrame(Delta = [result["Delta"]], Cpv = [result["Cr"]], Cadd_rate = result["Cadd_rate"], effectiveness = [result["Effectiveness"]], FES = [result["FES"]], Fim = [result["Fim"]], Fpv = [result["Fpv"]], SAP = [result["SAP"]]) |> CSV.write("imperfect_result$(episode).csv", append=true)
                    end
                end
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
