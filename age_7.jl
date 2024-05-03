function age_7_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    solve the household problem in period of child investment at age 7
    """
    @unpack c_size, a_size, h_size, h_grid, s_size = grids

    choices_initial = [0.5, 0.5, 0.5]
    choices_lb = [0.0, 0.0, 0.0]
    choices_ub = [1.0, 1.0, 1.0]

    # V_7 = zeros(h_size, a_size, s_size, h_size, a_size, c_size)
    for c_i = 1:c_size, a_i = 1:a_size, h_i = 1:h_size, s_i = 1:s_size, a_k_i = 1:a_size, h_k_i = 1:h_size
        states = [h_k_i, a_k_i, s_i, h_i, a_i, c_i]
        res = optimize(choices -> age_7_obj_function(choices, states, variables, prices, parameters, grids), choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()))
        @inbounds variables.V_7[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = -res.minimum
        @inbounds variables.policy_n_7[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[1]
        @inbounds variables.policy_l_2[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[2]
        @inbounds variables.policy_s_8[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[3]
        # println(states)
    end
end

function age_7_obj_function(choices::Vector{Float64}, states::Vector{Int64}, variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    construct the objective function at age 7
        choices: (n_7, l_2, s_8)
        states: (h_k_i, a_k_i, s_i, h_i, a_i, c_i)
    """
    @unpack β, b, q, γ_2, d_2, ω_2, ζ = parameters
    @unpack h_grid, a_grid, s_grid, s_min = grids

    # check if the number of choices is correct
    if length(choices) != 3
        error("dimension mismatch!")
    end

    # compute implied value for the given choices
    h_k_i, a_k_i, s_i, h_i, a_i, c_i = states
    n_7 = choices[1]
    l_2 = choices[2]
    l_2 = l_2 * (1.0 - n_7)
    h_7 = h_grid[h_i]
    a_7 = a_grid[a_i]
    earnings_7 = prices.w_S[c_i] * h_7 * (1.0 - n_7 - l_2)
    h_8 = a_7 * (n_7 * h_7)^b + h_7
    s_7 = s_grid[s_i]
    m_2 = prices.w_S[c_i] * h_7 * l_2 * (1.0 - γ_2) / γ_2
    x_2 = (γ_2 / prices.w_S[c_i])^γ_2 * (1.0 - γ_2)^(1.0 - γ_2) * (prices.w_S[c_i] * h_7 * l_2 + m_2 + d_2)
    x_2 = ζ * x_2
    h_k_2 = h_grid[h_k_i]
    a_k_2 = a_grid[a_k_i]
    h_k_3 = x_2^ω_2 * h_k_2^(1.0 - ω_2)
    budget = max(0.0, f_function(earnings_7 - m_2, s_7, 7, parameters, prices) - s_min)
    s_8 = choices[3] * budget
    c_7 = budget - s_8
    s_8 = s_8 + s_min
    itp_V_8 = V_8_itp_function(s_8, h_8, h_k_3, states, variables, grids)
    return -q * utility_function(c_7, parameters) - β * itp_V_8
end

function V_8_itp_function(s_8::Real, h_8::Real, h_k_3::Real, states::Vector{Int64}, variables::Mutable_Variables, grids::NamedTuple)
    """
    interpolate value function at age 8
        s_8: savings made at age 7
        h_8: imputed human capital at age 7 prior to market luck shock
        h_k_3: imputed kid's human capital at age 7 prior to market luck shock
    """
    # V_8 = zeros(h_size, a_size, s_size, h_size, a_size, c_size)
    @unpack ϵ_size, ϵ_grid = grids
    h_k_i, a_k_i, s_i, h_i, a_i, c_i = states
    V_8_temp = zeros(ϵ_size, 3)
    V_8_h_k_temp = zeros(2)
    s_ind, s_wgt = locate_s_function(s_8, grids)
    h_k_ind, h_k_wgt = locate_h_kid_function(h_k_3, grids)
    for ϵ_i = 1:ϵ_size
        h_8_ϵ = ϵ_grid[ϵ_i, c_i] * h_8
        h_ind, h_wgt = locate_h_function(h_8_ϵ, grids)
        V_8_temp[ϵ_i, 1] = sum(s_wgt .* variables.V_8[h_k_ind[1], a_k_i, s_ind[1]:s_ind[2], h_ind[1], a_i, c_i])
        V_8_temp[ϵ_i, 2] = sum(s_wgt .* variables.V_8[h_k_ind[1], a_k_i, s_ind[1]:s_ind[2], h_ind[2], a_i, c_i])
        V_8_temp[ϵ_i, 1] = sum(h_wgt .* V_8_temp[ϵ_i, 1:2])
        V_8_temp[ϵ_i, 2] = sum(s_wgt .* variables.V_8[h_k_ind[2], a_k_i, s_ind[1]:s_ind[2], h_ind[1], a_i, c_i])
        V_8_temp[ϵ_i, 3] = sum(s_wgt .* variables.V_8[h_k_ind[2], a_k_i, s_ind[1]:s_ind[2], h_ind[2], a_i, c_i])
        V_8_temp[ϵ_i, 2] = sum(h_wgt .* V_8_temp[ϵ_i, 2:3])
    end
    V_8_h_k_temp[1] = sum(V_8_temp[:, 1]) / ϵ_size
    V_8_h_k_temp[2] = sum(V_8_temp[:, 2]) / ϵ_size
    return sum(h_k_wgt .* V_8_h_k_temp[:])
end