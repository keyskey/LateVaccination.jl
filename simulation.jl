include("epidemics.jl")

module Simulation
    using ..Society
    using ..Epidemics
    using Random
    using CSV
    using DataFrames
    using Statistics
    using Printf

    # Calculation for a pair of Cr and Effectiveness. Continue until fv reaches to equilibrium
    function one_season(society::SocietyType, beta::Float64, gamma::Float64, m::Float64 , cr::Float64, effectiveness::Float64, num_initial_i::Int, num_initial_v::Int)
        society = Epidemics.initialize_society(society, beta, gamma, m, cr, effectiveness, num_initial_i, num_initial_v)
        fs0, fim0, fi0, fr0, fv0, fnv0 = Society.count_fraction(society)
        days = society.elapse_days
        @printf("Cr: %.1f e: %.1f Step: 0 Days: %.2f Fs: %.4f Fim: %.4f Fi: %.4f Fr: %.4f Fv: %.4f Fnv: %.4f \n", cr, effectiveness, days, fs0, fim0, fi0, fr0, fv0, fnv0)

        # SIR dynamics loop
        for timestep = 1:10000000
            society = Epidemics.gillespie_timestep(society, beta, gamma, m, cr, effectiveness, timestep)
            days = Epidemics.count_elapse_days(society)
            global fs, fim, fi, fr, fv, fnv = Society.count_fraction(society)
            fsv = society.num_sv/society.total_population
            @printf("Cr: %.1f e: %.1f Step: %i Days: %.2f Fs: %.4f Fim: %.4f Fi: %.4f Fr: %.4f Fv: %.4f Fnv: %.4f \n", cr, effectiveness, timestep, days, fs, fim, fi, fr, fv, fnv)

            # Check conversion
            if society.num_i == 0
                break
            end
        end

        SAP = Society.count_SAP(society)
        fr, fv, SAP
    end

    # Get data for one Cr-Effectiveness phase diagram
    function one_episode(total_population::Int, episode::Int)
        Random.seed!()
        DataFrame(Cr = [], Effectiveness = [], FES = [], VC = [], SAP = []) |> CSV.write("result$episode.csv")
        beta::Float64 = 0.000086  # 0.00008599
        gamma::Float64 = 1/3
        m::Float64 = 0.6
        num_initial_i::Int = 5
        num_initial_v::Int = div(total_population, 2)

        for cr::Float64 in [0.0] # = 0:0.1:1
            for effectiveness::Float64 in [1.0] # = 0:0.1:1
                society = SocietyType(total_population)
                FES, VC, SAP = one_season(society, beta, gamma, m, cr, effectiveness,  num_initial_i, num_initial_v)
                DataFrame(Cr = [cr], Effectiveness = [effectiveness], FES = [FES], VC = [VC], SAP = [SAP]) |> CSV.write("result$episode.csv", append=true)
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
