function age_9_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    solve the financial transfer problem of households
    """
    @unpack ϵ_size, ϵ_grid, c_size, a_size, a_grid, h_size, h_grid, s_size, s_grid, s_min, s_k_min, s_k_grid = grids
    @unpack β, θ, b = parameters

    choices_candidate = collect(0.1:0.2:0.9)
    states = zeros(5) # h_9, a_9, s_9, w_c, ϵ_c_mean
    @inbounds @views V_10_itp = linear_interpolation((s_grid, h_grid), zeros(s_size, h_size), extrapolation_bc=Interpolations.Flat())
    @inbounds @views V_4_itp = linear_interpolation(s_k_grid, zeros(s_size), extrapolation_bc=Interpolations.Flat())

    choices_initial = [0.5, 0.5, 0.5] # n_9, s_10, s_4
    choices_current = [0.5, 0.5, 0.5]
    choices_lb = [0.0, 0.0, 0.0]
    choices_ub = [1.0, 1.0, 1.0]
    function age_9_obj_function(choices::Vector{Float64}, states::Vector{Float64})
        """
        construct the objective function in age 9
            choices: (n_9, s_10, s_4)
            states: (h_k_4, a_k, c_k, s_9, h_9, a, c)
        """
        h_9, a_9, s_9, w_c, ϵ_c_mean = states
        n_9 = choices[1]
        h_10 = a_9 * (n_9 * h_9)^b + h_9
        earnings_9 = w_c * h_9 * (1.0 - n_9)
        budget = max(0.0, f_function(earnings_9, s_9, 9, parameters, prices) - s_min - s_k_min)
        s_10 = choices[2] * budget
        s_4 = choices[3] * (budget - s_10)
        c_9 = budget - s_10 - s_4
        s_10 = s_10 + s_min
        s_4 = s_4 + s_k_min
        return -utility_function(c_9, parameters) - θ * V_4_itp(s_4) - β * V_10_itp(s_10, ϵ_c_mean * h_10)
    end
    f = OptimizationFunction(age_9_obj_function)

    for c_i = 1:c_size
        @inbounds w_c = prices.w_S[c_i]
        @inbounds @views ϵ_c = ϵ_grid[:, c_i]
        ϵ_c_mean = mean(ϵ_c)
        @inbounds @views V_10_itp.itp.coefs .= variables.V_10[:, :, c_i]

        for c_k_i = 1:c_size, a_k_i = 1:a_size, h_k_i = 1:h_size
            @inbounds @views V_4_itp.itp.coefs .= variables.V_4[:, h_k_i, a_k_i, c_k_i]

            for a_i = 1:a_size, h_i = 1:h_size, s_i = 1:s_size
                # println("c_k_i = $c_k_i, a_k_i = $a_k_i, h_k_i = $h_k_i, a_i = $a_i, h_i = $h_i, s_i = $s_i")
                @inbounds h_9 = h_grid[h_i]
                @inbounds a_9 = a_grid[a_i]
                @inbounds s_9 = s_grid[s_i]
                states .= [h_9, a_9, s_9, w_c, ϵ_c_mean]
                for ch_1_i in choices_candidate, ch_2_i in choices_candidate, ch_3_i in choices_candidate
                    choices_current .= [ch_1_i, ch_2_i, ch_3_i]
                    if age_9_obj_function(choices_current, states) < age_9_obj_function(choices_initial, states)
                    choices_initial .= choices_current
                    end
                end
                prob = Optimization.OptimizationProblem(f, choices_initial, states, lb = choices_lb, ub = choices_ub)
                sol = solve(prob, NLopt.LN_SBPLX(), maxiters = 50, maxtime = 0.5)
                @inbounds variables.V_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = -sol.objective
                @inbounds variables.policy_n_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = sol.u[1]
                @inbounds variables.policy_s_10[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = sol.u[2]
                @inbounds variables.policy_s_4[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = sol.u[3]
            end
        end
    end
    return nothing
end