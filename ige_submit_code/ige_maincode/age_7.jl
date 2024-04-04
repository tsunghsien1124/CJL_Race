function age_7_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple)
    """
    solve the household problem in period of child investment
    """
    @unpack c_size, a_size, h_size, h_grid, s_size = parameters

    # V_7 = zeros(h_size, a_size, s_size, h_size, a_size, c_size)
    for c_i = 1:c_size, a_i = 1:a_size, h_i = 1:h_size, s_i = 1:s_size, a_k_i = 1:a_size, h_k_i = 1:h_size
        states = [h_k_i, a_k_i, s_i, h_i, a_i, c_i]
        choices_initial = [0.5, 0.5, 0.5]
        choices_lb = [0.0, 0.0, 0.0]
        choices_ub = [1.0, 1.0, 1.0]
        res = optimize(choices -> age_8_obj_function(choices, states, variables, prices, parameters), choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()))
        @inbounds variables.W_8[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = -res.minimum
        @inbounds variables.policy_n_8[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[1]
        @inbounds variables.policy_s_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[2]
        @inbounds variables.policy_n_3[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[3]
        # println(states)
    end
end

function age_7_obj_function(choices::Vector{Float64}, states::Vector{Int64}, variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple)
    """
    construct the objective function at age 7
        choices: (n_7, s_8, l_2)
        states: (h_k_i, a_k_i, s_i, h_i, a_i, c_i)
    """
    @unpack h_grid, a_grid, s_grid, β, b, s_min, q, γ_2, d_2 = parameters

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
    earnings_7 = prices.w_S[c_i] * h_7
    h_8 = a_7 * (n_7 * h_7)^b + h_7
    s_7 = s_grid[s_i]

    m_2 = earnings_7 * (1.0-γ_2) / γ_2 
    x_2 = (γ_2 / prices.w_S[c_i])^γ_2 * (1.0 - γ_2)^(1.0 - γ_2) * (earnings_7*l_2 + m_2 + d_2)

    h_k_3 = h_grid[h_k_i]
    a_k_3 = a_grid[a_k_i]

    mx = wp*lx *(1d0-kamk)/kamk
	xx = lamk*(wp*lx + mx + ddk)
	xx = zeta* xx

	hkx = xx**omegk * hk**(1d0-omegk)

    budget = max(0.0, f_function(earnings_8, s_8, 8, parameters) + f_function(earnings_k_3, 0.0, 3, parameters) - s_min)
    s_9 = choices[2] * budget
    c_8 = budget - s_9
    s_9 = s_9 + s_min
    itp_V_9 = V_9_itp_function(s_9, h_9, h_k_4, states, variables, parameters)
    return -q * utility_function(c_8, parameters) - β * itp_V_9
end

function V_9_itp_function(s_9::Real, h_9::Real, h_k_4::Real, states::Vector{Int64}, variables::Mutable_Variables, parameters::NamedTuple)
    """
    interpolate value function at age 9
        s_9: savings made at age 8
        h_9: imputed human capital at age 8 prior to market luck shock
        h_k_4: imputed kid's human capital at age 8 prior to market luck shock
    """
    # V_9 = zeros(h_size, a_size, c_size, s_size, h_size, a_size, c_size)
    @unpack ϵ_size, ϵ_grid = parameters
    h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i = states
    V_9_h_temp = zeros(ϵ_size, 3)
    V_9_h_k_temp = zeros(ϵ_size, 2)
    s_ind, s_wgt = locate_s_function(s_9, parameters)
    for ϵ_k_i = 1:ϵ_size
        h_k_4_ϵ = ϵ_grid[ϵ_k_i, c_k_i] * h_k_4
        h_k_ind, h_k_wgt = locate_h_function(h_k_4_ϵ, parameters)
        for ϵ_i = 1:ϵ_size
            h_9_ϵ = ϵ_grid[ϵ_i, c_i] * h_9
            h_ind, h_wgt = locate_h_function(h_9_ϵ, parameters)
            V_9_h_temp[ϵ_i, 1] = sum(s_wgt .* variables.V_9[h_k_ind[1], a_k_i, c_k_i, s_ind[1]:s_ind[2], h_ind[1], a_i, c_i])
            V_9_h_temp[ϵ_i, 2] = sum(s_wgt .* variables.V_9[h_k_ind[1], a_k_i, c_k_i, s_ind[1]:s_ind[2], h_ind[2], a_i, c_i])
            V_9_h_temp[ϵ_i, 1] = sum(h_wgt .* V_9_h_temp[ϵ_i, 1:2])
            V_9_h_temp[ϵ_i, 2] = sum(s_wgt .* variables.V_9[h_k_ind[2], a_k_i, c_k_i, s_ind[1]:s_ind[2], h_ind[1], a_i, c_i])
            V_9_h_temp[ϵ_i, 3] = sum(s_wgt .* variables.V_9[h_k_ind[2], a_k_i, c_k_i, s_ind[1]:s_ind[2], h_ind[2], a_i, c_i])
            V_9_h_temp[ϵ_i, 2] = sum(h_wgt .* V_9_h_temp[ϵ_i, 2:3])
        end
        V_9_h_k_temp[ϵ_k_i, 1] = sum(V_9_h_temp[:, 1]) / ϵ_size
        V_9_h_k_temp[ϵ_k_i, 2] = sum(V_9_h_temp[:, 2]) / ϵ_size
        V_9_h_k_temp[ϵ_k_i, 1] = sum(h_k_wgt .* V_9_h_k_temp[ϵ_k_i, :])
    end
    return sum(V_9_h_k_temp[:, 1]) / ϵ_size
end