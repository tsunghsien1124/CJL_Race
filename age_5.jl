function age_5_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    solve the household problem in period of child investment at age 5
    """
    @unpack c_size, a_size, a_grid, h_size, h_grid, h_k_grid, s_size, s_grid, s_min, ϵ_size, ϵ_grid = grids
    @unpack β, b, q, γ_0, d_0, ζ = parameters

    choices_candidate = collect(0.1:0.2:0.9)
    states = zeros(5) # h_5, a_5, s_5, w_c, ϵ_c_mean
    @inbounds @views V_6_itp = linear_interpolation((h_k_grid, s_grid, h_grid), zeros(h_size, s_size, h_size), extrapolation_bc=Interpolations.Flat())

    choices_initial = [0.5, 0.5, 0.5] # n_5, l_0, s_6
    choices_current = [0.5, 0.5, 0.5]
    choices_lb = [0.0, 0.0, 0.0]
    choices_ub = [1.0, 1.0, 1.0]
    function age_5_obj_function(choices::Vector{Float64}, states::Vector{Float64})
        """
        construct the objective function at age 5
            choices: (n_5, l_0, s_6)
            states: (s_5, h_5, a, c)
        """
        h_5, a_5, s_5, w_c, ϵ_c_mean = states
        n_5 = choices[1]
        l_0 = choices[2]
        l_0 = l_0 * (1.0 - n_5)
        earnings_5 = w_c * h_5 * (1.0 - n_5 - l_0)
        h_6 = a_5 * (n_5 * h_5)^b + h_5
        m_0 = w_c * h_5 * l_0 * (1.0 - γ_0) / γ_0
        x_0 = (γ_0 / w_c)^γ_0 * (1.0 - γ_0)^(1.0 - γ_0) * (w_c * h_5 * l_0 + m_0 + d_0)
        h_k_1 = ζ * x_0
        budget = max(0.0, f_function(earnings_5 - m_0, s_5, 5, parameters, prices) - s_min)
        s_6 = choices[3] * budget
        c_5 = budget - s_6
        s_6 = s_6 + s_min
        return -q * utility_function(c_5, parameters) - β * V_6_itp(h_k_1, s_6, ϵ_c_mean * h_6)
    end
    f = OptimizationFunction(age_5_obj_function)

    for c_i = 1:c_size
        @inbounds w_c = prices.w_S[c_i]
        @inbounds @views ϵ_c = ϵ_grid[:, c_i]
        ϵ_c_mean = mean(ϵ_c)
        for a_i = 1:a_size, a_k_i = 1:a_size
            @inbounds a_5 = a_grid[a_i]
            @inbounds @views V_6_itp.itp.coefs .= variables.V_6[:, a_k_i, :, :, a_i, c_i]
            for h_i = 1:h_size, s_i = 1:s_size
                @inbounds h_5 = h_grid[h_i]
                @inbounds s_5 = s_grid[s_i]
                states .= [h_5, a_5, s_5, w_c, ϵ_c_mean]
                for ch_1_i in choices_candidate, ch_2_i in choices_candidate, ch_3_i in choices_candidate
                    choices_current .= [ch_1_i, ch_2_i, ch_3_i]
                    if age_5_obj_function(choices_current, states) < age_5_obj_function(choices_initial, states)
                    choices_initial .= choices_current
                    end
                end
                prob = Optimization.OptimizationProblem(f, choices_initial, states, lb = choices_lb, ub = choices_ub)
                sol = solve(prob, NLopt.LN_SBPLX(), maxiters = 50, maxtime = 0.5)                
                @inbounds variables.V_5[a_k_i, s_i, h_i, a_i, c_i] = -sol.objective
                @inbounds variables.policy_n_5[a_k_i, s_i, h_i, a_i, c_i] = sol.u[1]
                @inbounds variables.policy_l_0[a_k_i, s_i, h_i, a_i, c_i] = sol.u[2]
                @inbounds variables.policy_s_6[a_k_i, s_i, h_i, a_i, c_i] = sol.u[3]
            end
        end
    end
    return nothing
end
