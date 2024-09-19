function age_9_function_th_1!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    solve the financial transfer problem of households
    """
    @unpack ϵ_size, ϵ_grid, c_size, a_size, a_grid, h_size, h_grid, s_size, s_grid, s_min, s_k_min, s_k_grid, loop_age_10 = grids
    @unpack β, θ, b = parameters

    choices_initial = [0.5, 0.5, 0.5]
    choices_lb = [0.0, 0.0, 0.0]
    choices_ub = [1.0, 1.0, 1.0]
    @inbounds @views V_10_itp = linear_interpolation((s_grid, h_grid), zeros(s_size, h_size), extrapolation_bc=Line())
    @inbounds @views V_4_itp = linear_interpolation(s_k_grid, zeros(s_size), extrapolation_bc=Line())

    for c_i = 1:c_size
        @inbounds w_c = prices.w_S[c_i]
        @inbounds @views ϵ_c = ϵ_grid[:, c_i]
        # @inbounds @views V_10_itp = linear_interpolation((s_grid, h_grid), variables.V_10[:, :, c_i], extrapolation_bc=Line())
        @inbounds @views V_10_itp.itp.coefs .= variables.V_10[:, :, c_i]
        @inbounds @views EV_10_itp(s_10, h_10) = (V_10_itp(s_10, ϵ_c[1] * h_10) + V_10_itp(s_10, ϵ_c[2] * h_10)) / ϵ_size
        for c_k_i = 1:c_size, a_k_i = 1:a_size, h_k_i = 1:h_size
            # @inbounds @views V_4_itp = linear_interpolation(s_k_grid, variables.V_4[:, h_k_i, a_k_i, c_k_i], extrapolation_bc=Line())
            @inbounds @views V_4_itp.itp.coefs .= variables.V_4[:, h_k_i, a_k_i, c_k_i]
            for a_i = 1:a_size, h_i = 1:h_size, s_i = 1:s_size
                @inbounds h_9 = h_grid[h_i]
                @inbounds a_9 = a_grid[a_i]
                @inbounds s_9 = s_grid[s_i]
                function age_9_obj_function(choices::Vector{Float64})
                    """
                    construct the objective function in age 9
                        choices: (n_9, s_10, s_4)
                        states: (h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i)
                    """
                    n_9 = choices[1]
                    h_10 = a_9 * (n_9 * h_9)^b + h_9
                    budget = max(0.0, f_function(w_c * h_9 * (1.0 - n_9), s_9, 9, parameters, prices) - s_min - s_k_min)
                    s_10 = choices[2] * budget
                    s_4 = choices[3] * (budget - s_10)
                    c_9 = budget - s_10 - s_4
                    s_10 = s_10 + s_min
                    s_4 = s_4 + s_k_min
                    return -utility_function(c_9, parameters) - θ * V_4_itp(s_4) - β * EV_10_itp(s_10, h_10)
                end
                res = optimize(age_9_obj_function, choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()))
                @inbounds variables.V_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = -res.minimum
                @inbounds variables.policy_n_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[1]
                @inbounds variables.policy_s_10[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[2]
                @inbounds variables.policy_s_4[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[3]
            end
        end
    end
    return nothing
end