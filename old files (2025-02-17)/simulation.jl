function simulation_function(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple; N::Int64=100, T::Int64=10)
    """
    simulate the economy with policy functions obtained with N hoseholds and T generations
    """
    @unpack a_size, a_grid, a_birth, a_Γ, age_size, ϵ_size, ϵ_grid, c_size, s_size, s_grid, s_k_grid, h_size, h_grid = grids
    @unpack s_min, s_max, s_k_min, s_k_max, h_min, h_max, h_k_max = grids
    @unpack b, γ_0, d_0, γ_1, d_1, γ_2, d_2, ζ, ω_1, ω_2, κ_tilde = parameters

    h_k_min = h_min

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
    s_simul[:, :, :] = rand(Categorical(ones(2) ./ 2), N * T * age_size)
    s_simul[:, 1, 1] = rand(Categorical(ones(s_size) ./ s_size), N)

    s_k_simul = zeros(Int, N, T)
    s_k_simul[:, :] = rand(Categorical(ones(2) ./ 2), N * T)

    h_simul = zeros(Int, N, T, age_size)
    h_simul[:, :, :] = rand(Categorical(ones(2) ./ 2), N * T * age_size)
    h_simul[:, 1, 1] = rand(Categorical(ones(h_size) ./ h_size), N)

    h_k_simul = zeros(Int, N, T, age_size)
    h_k_simul[:, :, :] = rand(Categorical(ones(2) ./ 2), N * T * age_size)

    for ind_i = 1:N
        for T_i = 1:(T-1)

            println("HH $ind_i in time $T_i")

            # (0) fixed over life cycle
            a_i = a_simul[ind_i, T_i]
            a_p_ = a_grid[a_i]
            c_i = c_simul[ind_i, T_i]
            w_p_ = prices.w_S[c_i]

            # (1) age 4
            h_4 = h_grid[h_simul[ind_i, T_i, 1]]
            s_4 = s_k_grid[s_simul[ind_i, T_i, 1]]
            w_p = w_p_ * h_4
            n_4 = variables.policy_n_4[s_simul[ind_i, T_i, 1], h_simul[ind_i, T_i, 1], a_i, c_i]
            n_4 = max(0.0, min(n_4, 1.0))
            earnings_4 = max(0.0, w_p * (1.0 - n_4))
            h_5 = ϵ_grid[ϵ_simul[ind_i, T_i, 2], c_i] * (a_p_ * (n_4 * h_4)^b + h_4)
            h_5 = max(h_min, min(h_5, h_max))
            budget = max(0.0, f_function(earnings_4, s_4, 4, parameters, prices) - s_min)
            s_5 = variables.policy_s_5[s_simul[ind_i, T_i, 1], h_simul[ind_i, T_i, 1], a_i, c_i] * budget
            s_5 = max(0.0, min(s_5, s_max)) + s_min

            # (2) age 5
            a_k_i = a_simul[ind_i, T_i+1]
            a_k_ = a_grid[a_k_i]

            s_ind, s_wgt = locate_s_function(s_5, grids)
            s_simul[ind_i, T_i, 2] = s_ind[s_simul[ind_i, T_i, 2]]
            s_5 = s_grid[s_simul[ind_i, T_i, 2]]

            h_ind, h_wgt = locate_h_function(h_5, grids)
            h_simul[ind_i, T_i, 2] = h_ind[h_simul[ind_i, T_i, 2]]
            h_5 = h_grid[h_simul[ind_i, T_i, 2]]
            w_p = w_p_ * h_5

            n_5 = variables.policy_n_5[a_k_i, s_simul[ind_i, T_i, 2], h_simul[ind_i, T_i, 2], a_i, c_i]
            l_0 = variables.policy_l_0[a_k_i, s_simul[ind_i, T_i, 2], h_simul[ind_i, T_i, 2], a_i, c_i]

            n_5 = max(0.0, min(n_5, 1.0))
            l_0 = l_0 * (1.0 - n_5)
            l_0 = max(0.0, min(l_0, 1.0 - n_5))
            earnings_5 = max(0.0, w_p * (1.0 - n_5 - l_0))

            h_6 = ϵ_grid[ϵ_simul[ind_i, T_i, 3], c_i] * (a_p_ * (n_5 * h_5)^b + h_5)
            h_6 = max(h_min, min(h_6, h_max))
            m_0 = w_p * l_0 * (1.0 - γ_0) / γ_0
            x_0 = (γ_0 / w_p_)^γ_0 * (1.0 - γ_0)^(1.0 - γ_0) * (w_p_ * h_5 * l_0 + m_0 + d_0)
            h_k_1 = ζ * x_0
            h_k_1 = max(h_k_min, min(h_k_1, h_k_max))

            budget = max(0.0, f_function(earnings_5 - m_0, s_5, 5, parameters, prices) - s_min)

            s_6 = variables.policy_s_6[a_k_i, s_simul[ind_i, T_i, 2], h_simul[ind_i, T_i, 2], a_i, c_i] * budget
            s_6 = max(0.0, min(s_6, s_max)) + s_min

            # (3) age 6
            s_ind, s_wgt = locate_s_function(s_6, grids)
            s_simul[ind_i, T_i, 3] = s_ind[s_simul[ind_i, T_i, 3]]
            s_6 = s_grid[s_simul[ind_i, T_i, 3]]

            h_ind, h_wgt = locate_h_function(h_6, grids)
            h_simul[ind_i, T_i, 3] = h_ind[h_simul[ind_i, T_i, 3]]
            h_6 = h_grid[h_simul[ind_i, T_i, 3]]
            w_p = w_p_ * h_6

            h_k_ind, h_k_wgt = locate_h_kid_function(h_k_1, grids)
            h_k_simul[ind_i, T_i, 3] = h_k_ind[h_k_simul[ind_i, T_i, 3]]
            h_k_1 = h_grid[h_k_simul[ind_i, T_i, 3]]

            n_6 = variables.policy_n_6[h_k_simul[ind_i, T_i, 3], a_k_i, s_simul[ind_i, T_i, 3], h_simul[ind_i, T_i, 3], a_i, c_i]
            l_1 = variables.policy_l_1[h_k_simul[ind_i, T_i, 3], a_k_i, s_simul[ind_i, T_i, 3], h_simul[ind_i, T_i, 3], a_i, c_i]
            l_1 = l_1 * (1.0 - n_6)
            earnings_6 = max(0.0, w_p * (1.0 - n_6 - l_1))
            h_7 = ϵ_grid[ϵ_simul[ind_i, T_i, 4], c_i] * (a_p_ * (n_6 * h_6)^b + h_6)
            h_7 = max(h_min, min(h_7, h_max))

            m_1 = w_p * l_1 * (1.0 - γ_1) / γ_1
            x_1 = (γ_1 / w_p_)^γ_1 * (1.0 - γ_1)^(1.0 - γ_1) * (w_p * l_1 + m_1 + d_1)
            x_1 = x_1 * ζ
            h_k_2 = x_1^ω_1 * h_k_1^(1.0 - ω_1)

            budget = max(0.0, f_function(earnings_6 - m_1, s_6, 6, parameters, prices) - s_min)

            s_7 = variables.policy_s_7[h_k_simul[ind_i, T_i, 3], a_k_i, s_simul[ind_i, T_i, 3], h_simul[ind_i, T_i, 3], a_i, c_i] * budget
            s_7 = max(0.0, min(s_7, s_max)) + s_min

            # (4) age 7
            s_ind, s_wgt = locate_s_function(s_7, grids)
            s_simul[ind_i, T_i, 4] = s_ind[s_simul[ind_i, T_i, 4]]
            s_7 = s_grid[s_simul[ind_i, T_i, 4]]

            h_ind, h_wgt = locate_h_function(h_7, grids)
            h_simul[ind_i, T_i, 4] = h_ind[h_simul[ind_i, T_i, 4]]
            h_7 = h_grid[h_simul[ind_i, T_i, 4]]
            w_p = w_p_ * h_7

            h_k_ind, h_k_wgt = locate_h_kid_function(h_k_2, grids)
            h_k_simul[ind_i, T_i, 4] = h_k_ind[h_k_simul[ind_i, T_i, 4]]
            h_k_2 = h_grid[h_k_simul[ind_i, T_i, 4]]

            n_7 = variables.policy_n_7[h_k_simul[ind_i, T_i, 4], a_k_i, s_simul[ind_i, T_i, 4], h_simul[ind_i, T_i, 4], a_i, c_i]
            l_2 = variables.policy_l_2[h_k_simul[ind_i, T_i, 4], a_k_i, s_simul[ind_i, T_i, 4], h_simul[ind_i, T_i, 4], a_i, c_i]
            l_2 = l_2 * (1.0 - n_7)
            earnings_7 = w_p * (1.0 - n_7 - l_2)
            h_8 = ϵ_grid[ϵ_simul[ind_i, T_i, 5], c_i] * (a_p_ * (n_7 * h_7)^b + h_7)
            h_8 = max(h_min, min(h_8, h_max))

            m_2 = w_p * l_2 * (1.0 - γ_2) / γ_2
            x_2 = (γ_2 / w_p_)^γ_2 * (1.0 - γ_2)^(1.0 - γ_2) * (w_p * l_2 + m_2 + d_2)
            x_2 = ζ * x_2
            h_k_3 = x_2^ω_2 * h_k_2^(1.0 - ω_2)

            budget = max(0.0, f_function(earnings_7 - m_2, s_7, 7, parameters, prices) - s_min)

            s_8 = variables.policy_s_8[h_k_simul[ind_i, T_i, 4], a_k_i, s_simul[ind_i, T_i, 4], h_simul[ind_i, T_i, 4], a_i, c_i] * budget
            s_8 = max(0.0, min(s_8, s_max)) + s_min

            # (5) age 8
            c_simul[ind_i, T_i+1] = variables.policy_c_3[h_k_simul[ind_i, T_i, 5], a_k_i, s_simul[ind_i, T_i, 5], h_simul[ind_i, T_i, 5], a_i, c_i] + 1
            c_k_i = c_simul[ind_i, T_i+1]
            w_k_ = prices.w_S[c_k_i]

            s_ind, s_wgt = locate_s_function(s_8, grids)
            s_simul[ind_i, T_i, 5] = s_ind[s_simul[ind_i, T_i, 5]]
            s_8 = s_grid[s_simul[ind_i, T_i, 5]]

            h_ind, h_wgt = locate_h_function(h_8, grids)
            h_simul[ind_i, T_i, 5] = h_ind[h_simul[ind_i, T_i, 5]]
            h_8 = h_grid[h_simul[ind_i, T_i, 5]]
            w_p = w_p_ * h_8

            h_k_ind, h_k_wgt = locate_h_kid_function(h_k_3, grids)
            h_k_simul[ind_i, T_i, 5] = h_k_ind[h_k_simul[ind_i, T_i, 5]]
            h_k_3 = h_grid[h_k_simul[ind_i, T_i, 5]]
            w_k = w_k_ * h_k_3

            n_8 = variables.policy_n_8[h_k_simul[ind_i, T_i, 5], a_k_i, c_k_i, s_simul[ind_i, T_i, 5], h_simul[ind_i, T_i, 5], a_i, c_i]

            earnings_8 = w_p * (1.0 - n_8)
            h_9 = ϵ_grid[ϵ_simul[ind_i, T_i, 6], c_i] * (a_p_ * (n_8 * h_8)^b + h_8)
            h_9 = max(h_min, min(h_9, h_max))

            n_k_3 = variables.policy_n_3[h_k_simul[ind_i, T_i, 5], a_k_i, c_k_i, s_simul[ind_i, T_i, 5], h_simul[ind_i, T_i, 5], a_i, c_i]
            n_k_3 = c_k_i == 2 ? κ_tilde + (1.0 - κ_tilde) * n_k_3 : n_k_3

            earnings_k_3 = w_k * (1.0 - n_k_3)

            h_k_4 = a_k_ * (n_k_3 * h_k_3)^b + h_k_3
            h_k_4 = ϵ_grid[ϵ_simul[ind_i, T_i+1, 1], c_k_i] * (a_k_ * (n_k_3 * h_k_3)^b + h_k_3)
            h_k_4 = max(h_k_min, min(h_k_4, h_k_max))

            budget = max(0.0, f_function(earnings_8, s_8, 8, parameters, prices) + f_function(earnings_k_3, 0.0, 3, parameters, prices) - s_min)

            s_9 = variables.policy_s_9[h_k_simul[ind_i, T_i, 5], a_k_i, c_k_i, s_simul[ind_i, T_i, 5], h_simul[ind_i, T_i, 5], a_i, c_i] * budget
            s_9 = max(0.0, min(s_9, s_max)) + s_min

            # (6) age 9
            s_ind, s_wgt = locate_s_function(s_9, grids)
            s_simul[ind_i, T_i, 6] = s_ind[s_simul[ind_i, T_i, 6]]
            s_9 = s_grid[s_simul[ind_i, T_i, 6]]

            h_ind, h_wgt = locate_h_function(h_9, grids)
            h_simul[ind_i, T_i, 6] = h_ind[h_simul[ind_i, T_i, 6]]
            h_9 = h_grid[h_simul[ind_i, T_i, 6]]
            w_p = w_p_ * h_9

            # h_k_ind, h_k_wgt = locate_h_kid_function(h_k_4, grids)
            # h_k_simul[ind_i, T_i, 6] = h_k_ind[h_k_simul[ind_i, T_i, 6]]
            # h_k_4 = h_grid[h_k_simul[ind_i, T_i, 6]]
            # w_k = w_k_ * h_k_4
            h_ind, h_wgt = locate_h_function(h_k_4, grids)
            h_simul[ind_i, T_i+1, 1] = h_ind[h_k_simul[ind_i, T_i, 6]]

            n_9 = variables.policy_n_9[h_k_simul[ind_i, T_i, 6], a_k_i, c_k_i, s_simul[ind_i, T_i, 6], h_simul[ind_i, T_i, 6], a_i, c_i]

            earnings_9 = w_p * (1.0 - n_9)
            h_10 = ϵ_grid[ϵ_simul[ind_i, T_i, 7], c_i] * (a_p_ * (n_9 * h_9)^b + h_9)
            h_10 = max(h_min, min(h_10, h_max))

            budget = max(0.0, f_function(earnings_9, s_9, 9, parameters, prices) - s_min - s_k_min)

            s_10 = variables.policy_s_10[h_k_simul[ind_i, T_i, 6], a_k_i, c_k_i, s_simul[ind_i, T_i, 6], h_simul[ind_i, T_i, 6], a_i, c_i]
            s_10 = s_10 * budget

            s_4 = variables.policy_s_4[h_k_simul[ind_i, T_i, 6], a_k_i, c_k_i, s_simul[ind_i, T_i, 6], h_simul[ind_i, T_i, 6], a_i, c_i]
            s_4 = s_4 * (budget - s_10)

            s_k_ind, s_k_wgt = locate_s_kid_function(s_4, grids)
            s_simul[ind_i, T_i+1, 1] = s_k_ind[s_k_simul[ind_i, T_i]]

            # (7) age 10
        end
    end

    return s_simul, h_simul
end