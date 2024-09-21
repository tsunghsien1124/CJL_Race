function age_6_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    solve the household problem in period of child investment at age 6
    """
    @unpack c_size, a_size, a_grid, h_size, h_grid, h_k_grid, s_size, s_grid, s_min, ϵ_size, ϵ_grid = grids
    @unpack β, b, q, γ_1, d_1, ω_1, ζ = parameters

    choices_initial = [0.5, 0.5, 0.5]
    choices_lb = [0.0, 0.0, 0.0]
    choices_ub = [1.0, 1.0, 1.0]
    @inbounds @views V_7_itp = linear_interpolation((h_k_grid, s_grid, h_grid), zeros(h_size, s_size, h_size), extrapolation_bc=Line())

    for c_i = 1:c_size
        @inbounds w_c = prices.w_S[c_i]
        @inbounds @views ϵ_c = ϵ_grid[:, c_i]
        ϵ_c_mean = mean(ϵ_c)
        for a_i = 1:a_size, a_k_i = 1:a_size
            @inbounds a_6 = a_grid[a_i]
            @inbounds a_k_1 = a_grid[a_k_i]
            @inbounds @views V_7_itp.itp.coefs .= variables.V_7[:, a_k_i, :, :, a_i, c_i]
            EV_7_itp(h_k_2, s_7, h_7) = V_7_itp(h_k_2, s_7, ϵ_c_mean * h_7)
            for h_i = 1:h_size, s_i = 1:s_size, h_k_i = 1:h_size
                @inbounds h_6 = h_grid[h_i]
                @inbounds s_6 = s_grid[s_i]
                @inbounds h_k_1 = h_k_grid[h_k_i]
                function age_6_obj_function(choices::Vector{Float64})
                    """
                    construct the objective function at age 6
                        choices: (n_6, l_1, s_7)
                        states: (h_k_i, a_k_i, s_i, h_i, a_i, c_i)
                    """
                    n_6 = choices[1]
                    l_1 = choices[2]
                    l_1 = l_1 * (1.0 - n_6)
                    earnings_6 = w_c * h_6 * (1.0 - n_6 - l_1)
                    h_7 = a_6 * (n_6 * h_6)^b + h_6
                    m_1 = w_c * h_6 * l_1 * (1.0 - γ_1) / γ_1
                    x_1 = (γ_1 / w_c)^γ_1 * (1.0 - γ_1)^(1.0 - γ_1) * (w_c * h_6 * l_1 + m_1 + d_1)
                    x_1 = x_1 * ζ
                    h_k_2 = x_1^ω_1 * h_k_1^(1.0 - ω_1)
                    budget = max(0.0, f_function(earnings_6 - m_1, s_6, 6, parameters, prices) - s_min)
                    s_7 = choices[3] * budget
                    c_6 = budget - s_7
                    s_7 = s_7 + s_min
                    return -q * utility_function(c_6, parameters) - β * EV_7_itp(h_k_2, s_7, h_7)
                end
                res = optimize(age_6_obj_function, choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()))
                @inbounds variables.V_6[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = -res.minimum
                @inbounds variables.policy_n_6[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[1]
                @inbounds variables.policy_l_1[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[2]
                @inbounds variables.policy_s_7[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[3]
            end
        end
    end
    return nothing
end