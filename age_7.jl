function age_7_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    solve the household problem in period of child investment at age 7
    """
    @unpack c_size, a_size, a_grid, h_size, h_grid, h_k_grid, s_size, s_grid, s_min, ϵ_size, ϵ_grid = grids
    @unpack β, b, q, γ_2, d_2, ω_2, ζ = parameters

    choices_initial = [0.5, 0.5, 0.5]
    choices_lb = [0.0, 0.0, 0.0]
    choices_ub = [1.0, 1.0, 1.0]
    @inbounds @views V_8_itp = linear_interpolation((h_k_grid, s_grid, h_grid), zeros(h_size, s_size, h_size), extrapolation_bc=Line())

    for c_i = 1:c_size
        @inbounds w_c = prices.w_S[c_i]
        @inbounds @views ϵ_c = ϵ_grid[:, c_i]
        ϵ_c_mean = mean(ϵ_c)
        for a_i = 1:a_size, a_k_i = 1:a_size
            @inbounds a_7 = a_grid[a_i]
            @inbounds a_k_2 = a_grid[a_k_i]
            @inbounds @views V_8_itp.itp.coefs .= variables.V_8[:, a_k_i, :, :, a_i, c_i]
            # @inbounds EV_8_itp(h_k_3, s_8, h_8) = (V_8_itp(h_k_3, s_8, ϵ_c[1] * h_8) + V_8_itp(h_k_3, s_8, ϵ_c[2] * h_8)) / ϵ_size
            EV_8_itp(h_k_3, s_8, h_8) = V_8_itp(h_k_3, s_8, ϵ_c_mean * h_8)
            for h_i = 1:h_size, s_i = 1:s_size, h_k_i = 1:h_size
                @inbounds h_7 = h_grid[h_i]
                @inbounds s_7 = s_grid[s_i]
                @inbounds h_k_2 = h_k_grid[h_k_i]
                function age_7_obj_function(choices::Vector{Float64})
                    """
                    construct the objective function at age 7
                        choices: (n_7, l_2, s_8)
                        states: (h_k_i, a_k_i, s_i, h_i, a_i, c_i)
                    """
                    n_7 = choices[1]
                    l_2 = choices[2]
                    l_2 = l_2 * (1.0 - n_7)
                    earnings_7 = w_c * h_7 * (1.0 - n_7 - l_2)
                    h_8 = a_7 * (n_7 * h_7)^b + h_7
                    m_2 = w_c * h_7 * l_2 * (1.0 - γ_2) / γ_2
                    x_2 = (γ_2 / w_c)^γ_2 * (1.0 - γ_2)^(1.0 - γ_2) * (w_c * h_7 * l_2 + m_2 + d_2)
                    x_2 = ζ * x_2
                    h_k_3 = x_2^ω_2 * h_k_2^(1.0 - ω_2)
                    budget = max(0.0, f_function(earnings_7 - m_2, s_7, 7, parameters, prices) - s_min)
                    s_8 = choices[3] * budget
                    c_7 = budget - s_8
                    s_8 = s_8 + s_min
                    return -q * utility_function(c_7, parameters) - β * EV_8_itp(h_k_3, s_8, h_8)
                end
                res = optimize(age_7_obj_function, choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()))
                @inbounds variables.V_7[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = -res.minimum
                @inbounds variables.policy_n_7[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[1]
                @inbounds variables.policy_l_2[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[2]
                @inbounds variables.policy_s_8[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = res.minimizer[3]
            end
        end
    end
    return nothing
end