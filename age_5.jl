function age_5_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple)
    """
    solve the household problem in period of child investment at age 5
    """
    @unpack c_size, a_size, h_size, h_grid, s_size = parameters

    # V_5 = zeros(a_size, s_size, h_size, a_size, c_size)
    for c_i = 1:c_size, a_i = 1:a_size, h_i = 1:h_size, s_i = 1:s_size, a_k_i = 1:a_size
        states = [a_k_i, s_i, h_i, a_i, c_i]
        choices_initial = [0.5, 0.5, 0.5]
        choices_lb = [0.0, 0.0, 0.0]
        choices_ub = [1.0, 1.0, 1.0]
        res = optimize(choices -> age_5_obj_function(choices, states, variables, prices, parameters), choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()))
        @inbounds variables.V_5[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = -res.minimum
        @inbounds variables.policy_n_5[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[1]
        @inbounds variables.policy_l_0[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[2]
        @inbounds variables.policy_s_6[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[3]
        # println(states)
    end
end

function age_5_obj_function(choices::Vector{Float64}, states::Vector{Int64}, variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple)
    """
    construct the objective function at age 5
        choices: (n_5, l_0, s_6)
        states: (a_k_i, s_i, h_i, a_i, c_i)
    """
    @unpack h_grid, a_grid, s_grid, β, b, s_min, q, γ_1, d_1, ω_1 = parameters

    # check if the number of choices is correct
    if length(choices) != 3
        error("dimension mismatch!")
    end

    # compute implied value for the given choices
    a_k_i, s_i, h_i, a_i, c_i = states
    n_6 = choices[1]
    l_1 = choices[2]
    l_1 = l_1 * (1.0 - n_6)
    h_6 = h_grid[h_i]
    a_6 = a_grid[a_i]
    earnings_6 = prices.w_S[c_i] * h_6 * (1.0 - n_6 - l_1)
    h_7 = a_6 * (n_6 * h_6)^b + h_6
    s_6 = s_grid[s_i]
    m_1 = prices.w_S[c_i] * h_6 * (1.0-γ_1) / γ_1 
    x_1 = (γ_1 / prices.w_S[c_i])^γ_1 * (1.0 - γ_1)^(1.0 - γ_1) * (prices.w_S[c_i] * h_6 *l_1 + m_1 + d_1)
    h_k_1 = h_grid[h_k_i]
    a_k_1 = a_grid[a_k_i]
    h_k_2 = x_1^ω_1 * h_k_1^(1.0-ω_1)
    budget = max(0.0, f_function(earnings_6-m_1, s_6, 6, parameters) - s_min)
    s_7 = choices[3] * budget
    c_6 = budget - s_7
    s_7 = s_7 + s_min
    itp_V_7 = V_7_itp_function(s_7, h_7, h_k_2, states, variables, parameters)
    return -q * utility_function(c_6, parameters) - β * itp_V_7
end

function V_7_itp_function(s_7::Real, h_7::Real, h_k_2::Real, states::Vector{Int64}, variables::Mutable_Variables, parameters::NamedTuple)
    """
    interpolate value function at age 7
        s_7: savings made at age 6
        h_7: imputed human capital at age 6 prior to market luck shock
        h_k_2: imputed kid's human capital at age 6 prior to market luck shock
    """
    # V_7 = zeros(h_size, a_size, s_size, h_size, a_size, c_size)
    @unpack ϵ_size, ϵ_grid = parameters
    h_k_i, a_k_i, s_i, h_i, a_i, c_i = states
    V_7_temp = zeros(ϵ_size, 3)
    V_7_h_k_temp = zeros(2)
    s_ind, s_wgt = locate_s_function(s_7, parameters)
    h_k_ind, h_k_wgt = locate_h_kid_function(h_k_2, parameters)
    for ϵ_i = 1:ϵ_size
        h_7_ϵ = ϵ_grid[ϵ_i, c_i] * h_7
        h_ind, h_wgt = locate_h_function(h_7_ϵ, parameters)
        V_7_temp[ϵ_i, 1] = sum(s_wgt .* variables.V_7[h_k_ind[1], a_k_i, s_ind[1]:s_ind[2], h_ind[1], a_i, c_i])
        V_7_temp[ϵ_i, 2] = sum(s_wgt .* variables.V_7[h_k_ind[1], a_k_i, s_ind[1]:s_ind[2], h_ind[2], a_i, c_i])
        V_7_temp[ϵ_i, 1] = sum(h_wgt .* V_7_temp[ϵ_i, 1:2])
        V_7_temp[ϵ_i, 2] = sum(s_wgt .* variables.V_7[h_k_ind[2], a_k_i, s_ind[1]:s_ind[2], h_ind[1], a_i, c_i])
        V_7_temp[ϵ_i, 3] = sum(s_wgt .* variables.V_7[h_k_ind[2], a_k_i, s_ind[1]:s_ind[2], h_ind[2], a_i, c_i])
        V_7_temp[ϵ_i, 2] = sum(h_wgt .* V_7_temp[ϵ_i, 2:3])
    end
    V_7_h_k_temp[1] = sum(V_7_temp[:, 1]) / ϵ_size
    V_7_h_k_temp[2] = sum(V_7_temp[:, 2]) / ϵ_size
    return sum(h_k_wgt .* V_7_h_k_temp[:])
end