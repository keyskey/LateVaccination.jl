module Decision
    using StatsBase
    using Match

    function choose_initial_pv(total_population)
        num_initial_pv::Int = div(total_population, 2)
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
    end

    function count_payoff(society, cr::Float64)
        for id in 1:society.total_population
            point::Float64 = @match society.state[id], society.strategy[id] begin
                "IM", _    => -cr
                "S" , "LV" => 0
                "R" , "LV" => -1
                 _  , _    => error("Error in count_payoff")               
            end
            society.point[id] = point
        end
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
    end
end
