include("society.jl")

module Epidemics
    using StatsBase
    using Random
    using Match
    using ..Society

    # Done
    function initialize_society(society::SocietyType, beta::Float64, gamma::Float64, m::Float64, cr::Float64, effectiveness::Float64, num_initial_i::Int, num_initial_v::Int)
        society = initialize_state(society, effectiveness, num_initial_i)
        society = initialize_strategy(society, num_initial_v)
        society = count_num_sv(society)
        society = reset_elapse_days(society) 
        society = set_initial_opponents(society)
        society = set_total_probability(society, beta, gamma, m, cr)
        
        society
    end

    # Done
    function initialize_state(society::SocietyType, effectiveness::Float64, num_initial_i::Int)
        society.state       = fill("S", society.total_population) 
        society.survivors   = Vector(1:society.total_population)  # Initially every agents are alive 
        society.point       = zeros(society.total_population)     # Need to count!!
        society.num_s       = society.total_population            # Determined
        society.num_im      = 0                                   # Determined
        society.num_i       = 0                                   # Need to count!!
        society.num_r       = 0                                   # Determined

        # このモデルでは最初からimmuned stateのエージェントは居ないので, initial iは誰から選んでも良い
        initial_i::Vector{Int} = StatsBase.self_avoid_sample!(Vector(1:society.total_population), Vector(1:num_initial_i))
        for i in initial_i
            society.state[i] = "I"
            society.num_s -= 1
            society.num_i += 1
            society.point[i] -= 1
        end

        society
    end
    
    # Done
    function initialize_strategy(society, num_initial_v::Int)
        society.num_v = 0
        society.num_nv = 0
        initial_v::Vector{Int} = StatsBase.self_avoid_sample!(Vector(1:society.total_population), Vector(1:num_initial_v))
        for i = 1:society.total_population
            if i in initial_v
                society.strategy[i] = "V"
                society.num_v += 1
            else
                society.strategy[i] = "NV"
                society.num_nv += 1
            end
        end

        society
    end
    
    # Done
    function count_num_sv(society)
        society.num_sv = 0
        for i = 1:society.total_population
            if society.state[i] == "S" && society.strategy[i] == "V"
                society.num_sv += 1
            end
        end

        society
    end
    
    # Done
    function reset_elapse_days(society::SocietyType)
        society.elapse_days = 0  # Determined

        society
    end

    # Done
    function set_initial_opponents(society::SocietyType)
        society.opponents = Random.shuffle(1:society.total_population)
        society
    end

    # Done
    function set_total_probability(society::SocietyType, beta::Float64, gamma::Float64, m::Float64, cr::Float64)
        # 状態遷移確率
        Ps_i = beta * society.num_i
        Pi_r = gamma

        # 戦略変更確率 (IBRA)
        # OPPはset_initial_opponentsで決定された初期対戦相手を使う(初期タイムステップ) or 
        # gillespie_timestepでstrategy_change loopの際に選んだ対戦相手をキャッシュしておいたものを使う(初期以降)
        P_strategy_change = 0
        for focal in society.survivors
            opp = society.opponents[focal]
            P_strategy_change += fermi(society.point[focal], society.point[opp])  # これの単位を/dayにすべく何かしらの係数をかける必要アリ
        end

        # ワクチン接種確率
        pi_sv_avg, pi_im_avg = averege_payoff(society)
        P_vaccinate = m * fermi(pi_sv_avg, pi_im_avg)

        # 全遷移確率
        society.total_probability = Ps_i * society.num_s + Pi_r * society.num_i + P_strategy_change + P_vaccinate * society.num_sv

        society
    end

    # Done
    function fermi(p_self::Float64, p_opp::Float64)
        1/(1 + exp( (p_self - p_opp)/0.1 ))
    end

    # Done V戦略側とNV戦略側の平均利得を返す
    # (=> SVグループとIMグループの平均利得比較に変更)
    # => やはりワクチン接種確率は感染者数に比例させないとCrへの感度が出ない
    # このコードではe < 1の時に一度ワクチンを打ったが免疫を獲得出来なかった人が再度ワクチンを打ちにいくことが許されている
    #    => ワクチンを打った段階でそのエージェントに"Vaccinated"フラグを付けといて、二度打ちに行かないようにする
    # やはり簡単化のためPerfect immunityモデルに絞りたい
    
    # Pending
    function gillespie_timestep(society::SocietyType, beta::Float64, gamma::Float64, m::Float64, cr::Float64, effectiveness::Float64, timestep::Int)
        rand_num::Float64 = rand()
        accum_probability::Float64 = 0.0
        someone_selected::Bool = false

        # State change loop(Done)
        for i in society.survivors
            accum_probability += @match society.state[i] begin
                "S" => beta * society.num_i/society.total_probability
                "I" => gamma/society.total_probability
                 _  => error("Error in state change, my state doesn't match with S/I")
            end

            if rand_num <= accum_probability
                society = state_change(society, i)
                someone_selected = true
                break  # Escape from survivors loop
            end
        end

        # Strategy change loop(IBRA)
        # total_probabilityの計算で使うOPPはここで決定してキャッシュしておく
        if someone_selected == false
            for focal in society.survivors
                if timestep == 1
                    opp = society.opponents[focal]
                    accum_probability += fermi(society.point[focal], society.point[opp])/society.total_probability
                    if rand_num <= accum_probability
                        society = strategy_change(society, focal, opp)
                        someone_selected = true
                        break  # Escape from survivors loop
                    end
                elseif timestep > 1
                    opp = rand(society.survivors)
                    if opp == focal
                        opp = rand(society.survivors)
                    end
                end

                society.opponents[focal] = opp  # Important: Cache opponent !!
                accum_probability += fermi(society.point[focal], society.point[opp])/society.total_probability
                if rand_num <= accum_probability
                    society = strategy_change(society, focal, opp)
                    someone_selected = true
                    break  # Escape from survivors loop
                end
            end
        end

        # Vaccination loop, only S state ∧ V strategy agents are selected
        if someone_selected == false
            for i in society.survivors
                if society.state[i] == "S" && society.strategy[i] == "V"
                    pi_sv_avg, pi_im_avg = averege_payoff(society)
                    accum_probability += m * fermi(pi_sv_avg, pi_im_avg)/society.total_probability
                    if rand_num <= accum_probability
                        society = vaccination(society, cr, effectiveness, i)
                        someone_selected = true
                        break  # Escape from survivors loop
                    end
                end
            end
        end

        # Error catch
        #=
        if someone_selected == false
            error("No one selected")
        end
        =#

        society = count_num_sv(society)
        society = set_total_probability(society, beta, gamma, m, cr)

        society
    end

    # Done
    function state_change(society::SocietyType, id::Int)
        @match society.state[id] begin
            "S" =>
              begin
                next_state = "I"
                society.num_s -= 1
                society.num_i += 1
                society.point[id] -= 1
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

        society
    end

    # Done
    function strategy_change(society::SocietyType, focal_id::Int, opp_id::Int)
        focal_strategy = society.strategy[focal_id]
        opp_strategy = society.strategy[opp_id]
        if focal_strategy == "V"
            if opp_strategy == "NV"
                society.num_v -= 1
                society.num_nv += 1
            end
        elseif focal_strategy == "NV" 
            if opp_strategy == "V"
                society.num_nv -= 1
                society.num_v += 1
            end
        end
        society.strategy[focal_id] = opp_strategy

        society
    end

    # Done
    function vaccination(society::SocietyType, cr::Float64, effectiveness::Float64, id::Int)
        println("Vaccination!")
        society.point[id] -= cr
        if rand() <= effectiveness
            society.state[id] = "IM"
            society.num_s -= 1
            society.num_im += 1
            society.survivors = filter(survivor_id -> survivor_id != id, society.survivors)  # If using filter!, survivors doesn't change outside of this begin-end block
        end

        society
    end

    # Done
    function count_elapse_days(society::SocietyType)
        society.elapse_days += log(1/rand())/society.total_probability

        society.elapse_days
    end
end
