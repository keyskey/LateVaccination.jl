# When running this file with multi-processing, type
# julia -p n main.jl  n: number of process
# after checking out the comment out at the beginning of each modules

@everywhere module Society
    using Statistics
    export SocietyType, Params

    mutable struct  SocietyType
        num_tot::Int
        num_s::Int
        num_v::Int
        num_im::Int
        num_i::Int
        num_r::Int
        state::Vector{Symbol}         # :S, :IM, :I, :R
        strategy::Vector{Symbol}      # :P or :L
        point::Vector{Float64}
        survivors::Vector{Int}
        total_probability::Float64
        elapse_days::Float64
        vaccinated::Vector{Bool}

        # Constructor
        SocietyType(population::Int) = new(
            population,                         # num_tot
            0,                                  # num_s
            0,                                  # num_v
            0,                                  # num_im
            0,                                  # num_i
            0,                                  # num_r
            fill(Symbol(""), population),       # state
            fill(Symbol(""), population),       # strategy
            zeros(population),                  # point
            Vector(1:population),               # survivors
            0,                                  # total_probability
            0,                                  # elapse_days
            fill(false, population)             # vaccinated
            )
    end

    struct Params
        beta::Float64
        gamma::Float64
        delta::Float64
        eps::Float64
        Cpv::Float64
        Clv::Float64
        Effectiveness::Float64
        num_initial_i::Int
        initial_pv::Vector{Int}
    end

    function count_state_fraction(society::SocietyType)
        fs = society.num_s/society.num_tot
        fv = society.num_v/society.num_tot
        fi = society.num_i/society.num_tot
        fr = society.num_r/society.num_tot

        return fs, fv, fi, fr
    end

    count_fpv(society::SocietyType) = length(society.strategy[society.strategy .== :P])/society.num_tot

    count_sap(society::SocietyType) = Statistics.mean(society.point)
end

@everywhere module Epidemics
    using StatsBase, Random, Printf
    using ..Society

    function initialize_state!(society::SocietyType, params::Params)
        society.survivors = Vector(1:society.num_tot)  # Need to remove vaccinated agents
        society.vaccinated = fill(false, society.num_tot)
        society.elapse_days = 0

        # Choose initial infected agent from LV
        LV_id = [id for (id, strategy) in enumerate(society.strategy) if strategy == :L]
        initial_i::Vector{Int} = StatsBase.self_avoid_sample!(LV_id, Vector(1:params.num_initial_i))

        society.num_s, society.num_v, society.num_im, society.num_i, society.num_r = 0, 0, 0, params.num_initial_i, 0

        for id in 1:society.num_tot
            if id in initial_i
                society.state[id] = :I
            elseif society.strategy[id] == :P
                society.vaccinated[id] = true
                society.num_v += 1
                if rand() < params.Effectiveness
                    society.state[id] = :IM
                    society.num_im += 1
                    deleteat!(society.survivors, findfirst(isequal(id), society.survivors))
                else
                    society.state[id] = :S
                    society.num_s += 1
                end
            else
                society.state[id] = :S
                society.num_s += 1
            end
        end

        set_total_probability!(society, params)
    end

    function set_total_probability!(society::SocietyType, params::Params)
        Ps_i = params.beta * society.num_i
        Pi_r = params.gamma
        P_v =  count_vaccination_probability(society, params)
        num_healthy_non_v = count_num_healthy_non_v(society)
        society.total_probability = Ps_i * society.num_s + Pi_r * society.num_i + P_v * num_healthy_non_v
    end

    count_num_healthy_non_v(society::SocietyType) = length([id for id in society.survivors if society.state[id] == :S && society.vaccinated[id] == false])

    count_vaccination_probability(society::SocietyType, params::Params) = params.delta * society.num_i / (params.Clv * society.num_v + params.eps)

    function one_season!(society::SocietyType, params::Params, season::Int)
        for timestep in 1:1000000
            rand_num = rand()
            accum_probability = 0.0
            go_to_state_change = true
            # someone_selected = false

            total_state_change_probability = params.beta * society.num_i * society.num_s + params.gamma * society.num_i
            if rand_num >= total_state_change_probability / society.total_probability
                accum_probability += total_state_change_probability
                go_to_state_change = false
            end

            # State change loop
            if go_to_state_change == true
                Ps_i = params.beta * society.num_i
                Pi_r = params.gamma
                for id in society.survivors
                    if society.state[id] == :S
                        accum_probability += Ps_i
                    elseif society.state[id] == :I
                        accum_probability += Pi_r
                    end
                    if rand_num <= accum_probability/society.total_probability
                        state_change!(society, id)
                        someone_selected = true
                        break
                    end
                end

            # Vaccination loop, only healthy non-vaccinators are selected
            else
                P_v = count_vaccination_probability(society, params)
                for id in society.survivors
                    if society.state[id] == :S && society.vaccinated[id] == false
                        accum_probability += P_v
                        if rand_num <= accum_probability/society.total_probability
                            vaccination!(society, id, params)
                            someone_selected = true
                            break
                        end
                    end
                end
            end

            set_total_probability!(society, params)
            days = count_elapse_days(society)

            """
            # Error catch
            if someone_selected == false
                error("Error in one_season, no one selected !!")
            end
            """

            # Check conversion
            if society.num_i == 0
                break
            end
        end
    end

    function state_change!(society::SocietyType, id::Int)
        if society.state[id] == :S
            next_state = :I
            society.num_s -= 1
            society.num_i += 1
        elseif society.state[id] == :I
            next_state =  :R
            society.num_i -= 1
            society.num_r += 1
            deleteat!(society.survivors, findfirst(isequal(id), society.survivors))
        end
        society.state[id] = next_state
    end

    function vaccination!(society::SocietyType, id::Int, params::Params)
        society.vaccinated[id] = true
        society.num_v += 1
        if rand() < params.Effectiveness
            society.state[id] = :IM
            society.num_s -= 1
            society.num_im += 1
            deleteat!(society.survivors, findfirst(isequal(id), society.survivors))
        end
    end

    function count_elapse_days(society::SocietyType)
        society.elapse_days += log(1/rand())/society.total_probability

        return society.elapse_days
    end
end

@everywhere module Decision
    using StatsBase
    using ..Society

    function choose_initial_pv(population::Int, num_initial_pv::Int)
        initial_pv::Vector{Int} = StatsBase.self_avoid_sample!(Vector(1:population), Vector(1:num_initial_pv))

        return initial_pv
    end

    function initialize_strategy!(society::SocietyType, params::Params)
        society.strategy = map(id -> ifelse(id in params.initial_pv, :P, :L), 1:society.num_tot)
    end

    function count_payoff!(society::SocietyType, params::Params)
        Cpv = params.Cpv
        Clv = params.Clv
        for id in 1:society.num_tot
            if society.strategy[id] == :P
                # Healthy Pre-emptive Vaccinators
                if society.state[id] in (:S, :IM)
                    society.point[id] = -Cpv
                # Infected Pre-emptive Vaccinators
                elseif society.state[id] == :R
                    society.point[id] = -Cpv-1
                end
            elseif society.strategy[id] == :L
                # Healthy Late Vaccinators
                if society.state[id] == :IM || (society.state[id] == :S && society.vaccinated[id] == true)
                    society.point[id] = -Clv
                # Healthy Non-Vaccinators
                elseif society.state[id] == :S && society.vaccinated[id] == false
                    society.point[id] = 0
                # Infected Late Vaccinators
                elseif society.state[id] == :R && society.vaccinated[id] == true
                    society.point[id] = -Clv-1
                # Infected non-vaccinators
                elseif society.state[id] == :R && society.vaccinated[id] == false
                    society.point[id] = -1
                end
            end
        end
    end

    function update_strategy!(society::SocietyType)
        next_strategy = copy(society.strategy)
        for id in 1:society.num_tot
            opp_id = rand(1:society.num_tot)
            while opp_id == id
                opp_id = rand(1:society.num_tot)
            end

            if society.strategy[opp_id] != society.strategy[id] && rand() < 1/(1 + exp((society.point[id] - society.point[opp_id])/0.1))
                next_strategy[id] = society.strategy[opp_id]
            end
        end

        society.strategy = copy(next_strategy)
    end
end

@everywhere module Simulation
    using Distributed, Random, CSV, DataFrames, Statistics, Printf
    using ..Epidemics, ..Decision, ..Society

    # Calculation for a pair of (Cpv, Clv, Effectiveness).
    function season_loop(society::SocietyType, params::Params)
        ###################### Initialization #####################
        Decision.initialize_strategy!(society, params)
        fpv0 = Society.count_fpv(society)
        fes_hist, fv_hist, sap_hist, fpv_hist = [], [], [], [fpv0]
        season = 0
        ###########################################################

        while true
            ######################## Initialization ########################
            season += 1
            Epidemics.initialize_state!(society, params)
            fs0, fv0, fi0, fr0 = Society.count_state_fraction(society)
            fpv_this_season = Society.count_fpv(society)
            @printf("Cpv: %.1f Clv: %.1f Effectiveness: %.1f Season: %i Step: %i Days: %.1f Fs: %.4f Fv: %.4f Fi: %.4f Fr: %.4f \n",
                    params.Cpv,
                    params.Clv,
                    params.Effectiveness,
                    season,
                    0,
                    society.elapse_days,
                    fs0,
                    fv0,
                    fi0,
                    fr0)
            ################################################################

            # One season
            Epidemics.one_season!(society, params, season)
            Decision.count_payoff!(society, params)
            Decision.update_strategy!(society)

            # Get final state
            fs, fv, fi, fes = Society.count_state_fraction(society)
            sap = Society.count_sap(society)
            fpv_next_season = Society.count_fpv(society)
            push!(fes_hist, fes)
            push!(fv_hist, fv)
            push!(sap_hist, sap)
            push!(fpv_hist, fpv_next_season)

            # Check conversion
            if (1 - fpv_next_season) * society.num_tot <= params.num_initial_i  # If number of LVs is less than num_initial_i, you can't choose initial infected agent!
                # Use the value at the previous season as a final solution
                global FES = fes_hist[end]
                global Fv = fv_hist[end]
                global SAP = sap_hist[end]
                global Fpv = fpv_hist[end-1]
                break

            elseif (season >= 100 && abs(Statistics.mean(fpv_hist[season-99:season]) - fpv_next_season) < 0.001) || season == 500
                global FES = Statistics.mean(fes_hist[season-99:season])
                global Fv  = Statistics.mean(fv_hist[season-99:season])
                global SAP = Statistics.mean(sap_hist[season-99:season])
                global Fpv = Statistics.mean(fpv_hist[season-99:season])
                break
            end
        end

        @printf("Cpv: %.1f Clv: %.1f Effectiveness: %.1f Finished with FES: %.4f Fv: %.4f Fpv: %.4f SAP: %.3f \n",
                params.Cpv,
                params.Clv,
                params.Effectiveness,
                FES,
                Fv,
                Fpv,
                SAP)

        result = Dict("Cpv"           => params.Cpv,
                      "Clv"           => params.Clv,
                      "Effectiveness" => params.Effectiveness,
                      "FES"           => FES,
                      "Fv"            => Fv,
                      "Fpv"           => Fpv,
                      "SAP"           => SAP)

        return result
    end

    # Get data for one Cr-Effectiveness phase diagram
    function one_episode(population::Int, episode::Int)
        Random.seed!()
        beta = 0.000086
        gamma = 1/3
        delta = 0.25
        eps = 0.1
        num_initial_i = 5
        num_initial_pv = div(population, 2)
        initial_pv = Decision.choose_initial_pv(population, num_initial_pv)
        society = SocietyType(population)

        result_file(delta) = "../../result/clv-dependent-δmodel/corrected_ver/result$(episode)_δ_$(delta).csv"

        DataFrame(Cpv = [], Clv = [], Effectiveness = [], FES = [], Fv = [], Fpv = [], SAP = []) |> CSV.write(result_file(delta))

        for Effectiveness in (0.1, 0.5, 0.8, 1.0)
            for Clv in 0:0.1:1.0
                results::Vector{Dict{String, Float64}} = Distributed.pmap(Cpv -> Simulation.season_loop(society, Params(beta, gamma, delta, eps, Cpv, Clv, Effectiveness, num_initial_i, initial_pv)), 0:0.1:1.0)
                # results::Vector{Dict{String, Float64}} = map(Cpv -> Simulation.season_loop(society, Params(beta, gamma, delta, eps, Cpv, Clv, Effectiveness, num_initial_i, initial_pv)), 0:0.1:1.0)

                for result in results
                    DataFrame(Cpv           = [result["Cpv"]],
                              Clv           = [result["Clv"]],
                              Effectiveness = [result["Effectiveness"]],
                              FES           = [result["FES"]],
                              Fv            = [result["Fv"]],
                              Fpv           = [result["Fpv"]],
                              SAP           = [result["SAP"]]) |> CSV.write(result_file(delta), append=true)
                end
            end
        end
    end
end

using .Simulation

const num_episode = 50
const population = 10000

for episode = 1:num_episode
    Simulation.one_episode(population, episode)
end
