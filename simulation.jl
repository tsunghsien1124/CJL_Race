function simulation_function(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple; N::Int64=1000, T::Int64=20)
    """
    simulate the economy with policy functions obtained with N hoseholds and T generations
    """
    @unpack a_size, a_birth, a_Γ, age_size, ϵ_size, c_size, s_size, h_size = grids

    # set seed
    Random.seed!(1124) # my birthday :)

    # draw ability
    a_simul = zeros(Int, N, T)
    a_simul[:, 1] = rand(Categorical(a_birth), N)
    for T_i = 2:T
        for ind_i = 1:N
            a_simul[ind_i, T_i] = rand(Categorical(a_Γ[a_simul[ind_i, T_i-1], :]))
        end
    end

    # draw luck shocks
    ϵ_simul = zeros(Int, N, T, age_size)
    ϵ_simul = rand(Categorical(ones(ϵ_size) ./ ϵ_size), (N, T, age_size))

    # draw initial states
    c_simul = zeros(Int, N, T)
    c_simul[:, 1] = rand(Categorical(ones(c_size) ./ c_size), N)
    s_simul = zeros(Int, N, T, age_size)
    s_simul[:, 1, 1] = rand(Categorical(ones(s_size) ./ s_size), N)
    h_simul = zeros(Int, N, T, age_size)
    h_simul[:, 1, 1] = rand(Categorical(ones(h_size) ./ h_size), N)

    for ind_i = 1:N
        for T_i = 1:(T-2)

            # (1) age 4
            a_p = a_grid[a_simul[ind_i, T_i]]
            h_4 = h_grid[h_simul[ind_i, T_i, 1]]
            s_4 = s_k_grid[s_simul[ind_i, T_i, 1]]
            w_p = prices.w_S[c_simul[:, 1]] * h_4
            n_4 = variables.policy_n_4[s_simul[ind_i, T_i, 1], h_simul[ind_i, T_i, 1], a_simul[ind_i, T_i], c_simul[ind_i, 1]]
            n_4 = max(0.0, min(n_4, 1.0))
            earnings_4 = max(0.0, w_p * (1.0 - n_4))
            h_5 = ϵ_grid[ϵ_simul[ind_i, T_i, 2], c_simul[ind_i, 1]] * (a_p * (n_4 * h_4)^b + h_4)
            h_5 = max(h_min, min(h_5, h_max))
            budget = max(0.0, f_function(earnings_4, s_4, 4, parameters, prices) - s_min)
            s_5 = variables.policy_s_5[s_simul[ind_i, T_i, 1], h_simul[ind_i, T_i, 1], a_simul[ind_i, T_i], c_simul[ind_i, 1]] * budget
            s_5 = max(0.0, min(s_5, s_max)) + s_min

            # (2) age 5

            # (3) age 6

            # (4) age 7

            # (5) age 8

            # (6) age 9

            # (7) age 10

        end
    end

end