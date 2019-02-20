include("epidemics.jl")
include("decision.jl")

module Simulation
    using ..Epidemics
    using ..Decision
    using ..Society
    using Random
    using CSV
    using DataFrames
    using Statistics
    using Printf

    # Calculation for a pair of Cr and Effectiveness. Continue until fv reaches to equilibrium
    function season_loop(society::SocietyType, beta::Float64, gamma::Float64, delta::Float64 , cr::Float64, num_initial_i::Int, initial_pv::Vector{Int})
        # Initialization
        society = Decision.initialize_strategy(society, initial_pv)
        fpv0 = Society.count_strategy_fraction(society)
        fpv_hist = [fpv0]

        # Season loop
        for season in 1:500
            # Initialization
            society = Epidemics.initialize_epidemics(society, beta, gamma, delta, num_initial_i)
            fs0, fim0, fi0, fr0 = Society.count_state_fraction(society)
            fpv = Society.count_strategy_fraction(society)
            @printf("Cr: %.1f Season: %i Step: 0 Days: %.2f Fs: %.4f Fim: %.4f Fpv: %.4f Fi: %.4f Fr: %.4f \n", cr, season, society.elapse_days, fs0, fim0, fpv, fi0, fr0)
            # DataFrame(fs = [], fim = [], fi = [], fr = [], IoverIM = []) |> CSV.write("time_series_delta_$(delta)_cr_$(cr)_season_$(season).csv")

            # One season
            society = Decision.update_strategy(Decision.count_payoff(Epidemics.one_season(society, beta, gamma, delta, cr, season), cr))
            global fs, Fim, fi, FES = Society.count_state_fraction(society)
            global SAP = Society.count_SAP(society)
            fpv_next = Society.count_strategy_fraction(society)
            push!(fpv_hist, fpv_next)
            @printf("Cr: %.1f Next season start with: Fim: %.4f \n", cr, fpv_next)

            # Check conversion
            if fpv_next * society.total_population >= society.total_population - num_initial_i  # Can't choose initial infected agent!
                global Fpv = fpv
                break
            elseif season >= 100 && Statistics.mean(fpv_hist[season-99:season]) - fpv_next < 0.001
                global Fpv = Statistics.mean(fpv_hist[season-98:season+1])
                break
            end
        end

        @printf("Cr: %.1f Finished with FES: %.4f Fim(VC): %.4f Fpv: %.4f SAP: %.3f \n", cr, FES, Fim, Fpv, SAP)

        return FES, Fim, Fpv, SAP
    end

    # Get data for one Cr-Effectiveness phase diagram
    function one_episode(total_population::Int, episode::Int)
        Random.seed!()
        beta::Float64 = 0.000086  # 0.00008599
        gamma::Float64 = 1/3
        num_initial_i::Int = 5
        num_initial_pv::Int = div(total_population, 10)
        initial_pv::Vector{Int} = Decision.choose_initial_pv(total_population, num_initial_pv)
        society = SocietyType(total_population)

        for delta::Float64 in [0.2, 0.5, 0.8, 1.0]
            DataFrame(Cr = [], FES = [], Fim = [], Fpv = [], SAP = []) |> CSV.write("result$(episode)_delta_$(delta).csv")
            for cr::Float64 in 0:0.1:1
                FES, fim, fpv, SAP = season_loop(society, beta, gamma, delta, cr, num_initial_i, initial_pv)
                DataFrame(Cr = [cr], FES = [FES], Fim = [fim], Fpv = [fpv], SAP = [SAP]) |> CSV.write("result$(episode)_delta_$(delta).csv", append=true)
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
