function age_4_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple)
    """
    solve the household problem in period of independece at age 4
    """
    @unpack c_size, a_size, h_size, h_grid, s_size = parameters

    # V_4 = zeros(s_size, h_size, a_size, c_size)
    for c_i = 1:c_size, a_i = 1:a_size, h_i = 1:h_size, s_i = 1:s_size
        states = [s_i, h_i, a_i, c_i]
        choices_initial = [0.5, 0.5, 0.5]
        choices_lb = [0.0, 0.0, 0.0]
        choices_ub = [1.0, 1.0, 1.0]
        res = optimize(choices -> age_4_obj_function(choices, states, variables, prices, parameters), choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()))
        @inbounds variables.V_4[s_i, h_i, a_i, c_i] = -res.minimum
        @inbounds variables.policy_n_4[s_i, h_i, a_i, c_i] = res.minimizer[1]
        @inbounds variables.policy_s_5[s_i, h_i, a_i, c_i] = res.minimizer[2]
        # println(states)
    end
end

function age_4_obj_function(choices::Vector{Float64}, states::Vector{Int64}, variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple)
    """
    construct the objective function at age 4
        choices: (n_4, s_5)
        states: (s_i, h_i, a_i, c_i)
    """
    @unpack h_grid, a_grid, s_grid, β, b, s_min = parameters

    # check if the number of choices is correct
    if length(choices) != 3
        error("dimension mismatch!")
    end

    # compute implied value for the given choices
    s_i, h_i, a_i, c_i = states
    n_4 = choices[1]
    h_4 = h_grid[h_i]
    a_4 = a_grid[a_i]
    earnings_4 = prices.w_S[c_i] * h_4 * (1.0 - n_4)
    h_5 = a_4 * (n_4 * h_4)^b + h_4
    s_4 = s_grid[s_i]
    budget = max(0.0, f_function(earnings_4, s_4, 4, parameters) - s_min)
    s_5 = choices[2] * budget
    c_4 = budget - s_5
    s_5 = s_5 + s_min
    itp_V_5 = V_5_itp_function(s_5, h_5, states, variables, parameters)
    return -utility_function(c_4, parameters) - β * itp_V_5
end

function V_5_itp_function(s_5::Real, h_5::Real, states::Vector{Int64}, variables::Mutable_Variables, parameters::NamedTuple)
    """
    interpolated value function at age 5
        s_5: savings made at age 4
        h_5: imputed human capital at age 4 prior to market luck shock
    """
    # V_5 = zeros(a_size, s_size, h_size, a_size, c_size)
    @unpack ϵ_size, ϵ_grid, a_Γ, a_size = parameters
    s_i, h_i, a_i, c_i = states
    V_5_temp = zeros(ϵ_size, 2)
    V_5_a_temp = zeros(a_size)
    s_ind, s_wgt = locate_s_function(s_5, parameters)
    for a_k_i = 1:a_size
        for ϵ_i = 1:ϵ_size
            h_5_ϵ = ϵ_grid[ϵ_i, c_i] * h_5
            h_ind, h_wgt = locate_h_function(h_5_ϵ, parameters)
            V_5_temp[ϵ_i, 1] = sum(s_wgt .* variables.V_5[a_k_i, s_ind[1]:s_ind[2], h_ind[1], a_i, c_i])
            V_5_temp[ϵ_i, 2] = sum(s_wgt .* variables.V_5[a_k_i, s_ind[1]:s_ind[2], h_ind[2], a_i, c_i])
            V_5_temp[ϵ_i, 1] = sum(h_wgt .* V_5_temp[ϵ_i, 1:2])
        end
        V_5_a_temp[a_i] = sum(V_5_temp[:, 1]) / ϵ_size
    end
    return sum(a_Γ[a_i,:] .* V_5_a_temp[:])
end