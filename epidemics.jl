include("society.jl")

module Epidemics
    using StatsBase
    using Match
    using Random
    using Printf
    using ..Society

    function initialize_epidemics(society::SocietyType, beta::Float64, gamma::Float64, m::Float64, num_initial_i::Int)
        initialize_state(society, num_initial_i)
        reset_elapse_days(society) 
        set_total_probability(society, beta, gamma, m)
        
    end
    
    function reset_elapse_days(society::SocietyType)
        society.elapse_days = 0
    end

    function initialize_state(society::SocietyType, num_initial_i::Int)
        society.survivors = Vector(1:society.total_population)  # Need to remove IM state

        # Pre-emptive Vは最初からIM stateなのでinitial iはLVから選ぶ
        LV_id = [id for (id, strategy) in enumerate(society.strategy) if strategy == "LV"]
        initial_i::Vector{Int} = StatsBase.self_avoid_sample!(LV_id, Vector(1:num_initial_i))
        society.num_s       = 0  # Need to count!!
        society.num_im      = 0  # Need to count!!
        society.num_i       = 0  # Need to count!!
        society.num_r       = 0  # Determined
        for id in 1:society.total_population
            if id in initial_i
                society.state[id] = "I"
                society.num_i += 1
            elseif society.strategy[id] == "PV"
                society.state[id] = "IM"
                society.survivors = filter(survivor_id -> survivor_id != id, society.survivors)  # If using filter!, survivors doesn't change outside of this begin-end block
                society.num_im += 1
            else
                society.state[id] = "S"
                society.num_s += 1
            end
        end
    end
    
    function set_total_probability(society::SocietyType, beta::Float64, gamma::Float64, m::Float64)
        Ps_i = beta * society.num_i
        Pi_r = gamma
        P_vaccinate =  m * society.num_i/society.num_im
        society.total_probability = (Ps_i + P_vaccinate) * society.num_s + Pi_r * society.num_i
    end

    function one_season(society::SocietyType, beta::Float64, gamma::Float64, m::Float64, cr::Float64)
        for timestep in 1:10000000
            rand_num::Float64 = rand()
            accum_probability::Float64 = 0.0
            someone_selected::Bool = false

            # State change loop
            for id in society.survivors
                accum_probability += @match society.state[id] begin
                    "S" => beta * society.num_i/society.total_probability
                    "I" => gamma/society.total_probability
                    _  => error("Error in state change, my state doesn't match with S/I")
                end

                if rand_num <= accum_probability
                    state_change(society, id)
                    someone_selected = true
                    break  # Escape from survivors loop
                end
            end

            # Vaccination loop, only health free riders are selected
            if someone_selected == false
                for id in society.survivors
                    if society.state[id] == "S" && society.strategy[id] == "LV"
                        accum_probability += m*society.num_i/society.num_im/society.total_probability
                        if rand_num <= accum_probability
                            vaccination(society, id)
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

            set_total_probability(society, beta, gamma, m)
            days = count_elapse_days(society)
            global fs, fim, fi, fr = Society.count_state_fraction(society)
            @printf("Cr: %.1f Step: %i Days: %.2f Fs: %.4f Fim: %.4f Fi: %.4f Fr: %.4f \n", cr, timestep, days, fs, fim, fi, fr)

            # Check conversion
            if society.num_i == 0
                break
            end
        end
        
        return fr
    end

    function state_change(society::SocietyType, id::Int)
        @match society.state[id] begin
            "S" =>
              begin
                next_state = "I"
                society.num_s -= 1
                society.num_i += 1
              end

            "I" =>
              begin
                next_state = "R"
                society.num_i -= 1
                society.num_r += 1
                society.survivors = filter(survivor_id -> survivor_id != id, society.survivors)  # If using filter!, survivors doesn't change outside of this begin-end block
              end
            _ => error("Error in desease spreading, R or the other state is picked up in gillespie method")
        end
        society.state[id] = next_state
    end

    function vaccination(society::SocietyType, id::Int)
        society.state[id] = "IM"
        society.num_s -= 1
        society.num_im += 1
        society.survivors = filter(survivor_id -> survivor_id != id, society.survivors)  # If using filter!, survivors doesn't change outside of this begin-end block
    end

    function count_elapse_days(society::SocietyType)
        society.elapse_days += log(1/rand())/society.total_probability
    end
end
