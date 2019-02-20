include("society.jl")

module Epidemics
    using StatsBase
    using Random
    using Printf
    using DataFrames
    using CSV
    using ..Society

    function initialize_epidemics(society::SocietyType, beta::Float64, gamma::Float64, delta::Float64, num_initial_i::Int)
        society = set_total_probability(reset_elapse_days(initialize_state(society, num_initial_i)), beta, gamma, delta)

        return society
    end
    
    function reset_elapse_days(society::SocietyType)
        society.elapse_days = 0

        return society
    end

    function initialize_state(society::SocietyType, num_initial_i::Int)
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
    
    function set_total_probability(society::SocietyType, beta::Float64, gamma::Float64, delta::Float64)
        Ps_i = beta * society.num_i
        Pi_r = gamma
        P_vaccinate =  delta * society.num_i/(society.num_im+1)
        society.total_probability = (Ps_i + P_vaccinate) * society.num_s + Pi_r * society.num_i

        return society
    end

    function one_season(society::SocietyType, beta::Float64, gamma::Float64, delta::Float64, cr::Float64, season::Int)
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

            # Vaccination loop, only health free riders are selected
            if someone_selected == false
                for id in society.survivors
                    if society.state[id] == "S" && society.strategy[id] == "LV"
                        accum_probability += delta * society.num_i/(society.num_im+1)/society.total_probability
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
            days = count_elapse_days(society)
            fs, fim, fi, fr = Society.count_state_fraction(society)
            @printf("Cr: %.1f Season: %i Step: %i Days: %.2f Fs: %.4f Fim: %.4f Fi: %.4f Fr: %.4f \n", cr, season, timestep, days, fs, fim, fi, fr)
            # DataFrame(fs = [fs], fim = [fim], fi = [fi], fr = [fr], IoverIM = [society.num_i/society.num_im]) |> CSV.write("time_series_delta_$(delta)_cr_$(cr)_season_$(season).csv", append=true)

            # Check conversion
            if society.num_i == 0
                break
            end
        end
        
        return society
    end

    function state_change(society::SocietyType, id::Int)
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

    function vaccination(society::SocietyType, id::Int)
        society.state[id] = "IM"
        society.num_s -= 1
        society.num_im += 1
        deleteat!(society.survivors, findfirst(isequal(id), society.survivors))

        return society
    end

    function count_elapse_days(society::SocietyType)
        society.elapse_days += log(1/rand())/society.total_probability

        return society.elapse_days
    end
end
