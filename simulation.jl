include("epidemics.jl")
include("decision.jl")

module Simulation
    using ..Society
    using ..Epidemics
    using ..Decision
    using Random
    using CSV
    using DataFrames
    using Statistics
    using Printf

    # Calculation for a pair of Cr and Effectiveness. Continue until fv reaches to equilibrium
    function season_loop(society::SocietyType, beta::Float64, gamma::Float64, m::Float64 , cr::Float64, num_initial_i::Int, initial_pv::Vector{Int})
        Decision.initialize_strategy(society, initial_pv)

        fpv0, flv0 = Society.count_strategy_fraction(society)
        fpv_hist = [fpv0]

        for season in 1:1000
            Epidemics.initialize_epidemics(society, beta, gamma, m, num_initial_i)
            fs0, fim0, fi0, fr0 = Society.count_state_fraction(society)
            days = society.elapse_days
            @printf("Cr: %.1f Step: 0 Days: %.2f Fs: %.4f Fim: %.4f Fi: %.4f Fr: %.4f Fpv: %.4f Flv: %.4f \n", cr, days, fs0, fim0, fi0, fr0, fpv0, flv0)

            global FES = Epidemics.one_season(society, beta, gamma, m, cr)
            Decision.count_payoff(society, cr)
            Decision.update_strategy(society)
            global fpv, flv = Society.count_strategy_fraction(society)
            push!(fpv_hist, fpv)

            # Check conversion
            if season >= 100 && Statistics.mean(fpv_hist[season-99:season]) - fpv < 0.001
                fpv = Statistics.mean(fpv_hist[season-98:season+1])
                break
            end
        end

        SAP = Society.count_SAP(society)
        
        return FES, fpv, SAP
    end

    # Get data for one Cr-Effectiveness phase diagram
    function one_episode(total_population::Int, episode::Int)
        Random.seed!()
        DataFrame(Cr = [], FES = [], VC = [], SAP = []) |> CSV.write("result$(episode).csv")
        beta::Float64 = 0.000086  # 0.00008599
        gamma::Float64 = 1/3
        m::Float64 = 0.6
        num_initial_i::Int = 5
        initial_pv::Vector{Int} = Decision.choose_initial_pv(total_population)
        society = SocietyType(total_population)

        for cr::Float64 in [0.1] # = 0:0.1:1
            FES, fpv, SAP = season_loop(society, beta, gamma, m, cr, num_initial_i, initial_pv)
            DataFrame(Cr = [cr], FES = [FES], VC = [VC], SAP = [SAP]) |> CSV.write("result$(episode).csv", append=true)
        end
    end
end

using .Simulation

const num_episode = 1
const total_population = 10000

for episode = 1:num_episode
    Simulation.one_episode(total_population, episode)
end
