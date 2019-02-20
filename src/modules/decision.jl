module Decision
    using StatsBase

    function choose_initial_pv(total_population::Int, num_initial_pv::Int)
        initial_pv::Vector{Int} = StatsBase.self_avoid_sample!(Vector(1:total_population), Vector(1:num_initial_pv))

        return initial_pv
    end

    function initialize_strategy(society, initial_pv::Vector{Int})
        for id in 1:society.total_population
            if id in initial_pv
                society.strategy[id] = "PV"
            else
                society.strategy[id] = "LV"
            end
        end

        return society
    end

    function count_payoff(society, cr::Float64)
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

    # IBRA
    function update_strategy(society)
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
