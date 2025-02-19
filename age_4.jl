function age_4_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    solve the household problem in period of independece at age 4
    """
    @unpack c_size, a_size, a_grid, a_Γ, h_size, h_grid, s_size, s_grid, s_min, ϵ_size, ϵ_grid = grids
    @unpack β, b = parameters

    choices_candidate = collect(0.1:0.2:0.9)
    states = zeros(5) # h_4, a_4, s_4, w_c, ϵ_c_mean
    @inbounds @views V_5_itp = linear_interpolation((s_grid, h_grid), zeros(s_size, h_size), extrapolation_bc=Line())
    EV_5 = zeros(s_size, h_size)

    choices_initial = [0.5, 0.5] # n_4, s_5
    choices_current = [0.5, 0.5]
    choices_lb = [0.0, 0.0]
    choices_ub = [1.0, 1.0]
    function age_4_obj_function(choices::Vector{Float64}, states::Vector{Float64})
        """
        construct the objective function at age 4
            choices: (n_4, s_5)
            states: (s_4, h_4, a, c)
        """
        h_4, a_4, s_4, w_c, ϵ_c_mean = states
        n_4 = choices[1]
        earnings_4 = w_c * h_4 * (1.0 - n_4)
        h_5 = a_4 * (n_4 * h_4)^b + h_4
        budget = max(0.0, f_function(earnings_4, s_4, 4, parameters, prices) - s_min)
        s_5 = choices[2] * budget
        c_4 = budget - s_5
        s_5 = s_5 + s_min
        return -utility_function(c_4, parameters) - β * V_5_itp(s_5, ϵ_c_mean * h_5)
    end
    f = OptimizationFunction(age_4_obj_function)

    for c_i = 1:c_size
        @inbounds w_c = prices.w_S[c_i]
        @inbounds @views ϵ_c = ϵ_grid[:, c_i]
        ϵ_c_mean = mean(ϵ_c)
        for a_i = 1:a_size
            @inbounds a_4 = a_grid[a_i]
            EV_5 .= 0.0
            for a_k_i = 1:a_size
                EV_5 .+= a_Γ[a_i, a_k_i] .* variables.V_5[a_k_i, :, :, a_i, c_i]
            end
            @inbounds @views V_5_itp.itp.coefs .= EV_5
            for h_i = 1:h_size, s_i = 1:s_size
                @inbounds h_4 = h_grid[h_i]
                @inbounds s_4 = s_grid[s_i]
                states .= [h_4, a_4, s_4, w_c, ϵ_c_mean]
                for ch_1_i in choices_candidate, ch_2_i in choices_candidate
                    choices_current .= [ch_1_i, ch_2_i]
                    if age_4_obj_function(choices_current, states) < age_4_obj_function(choices_initial, states)
                    choices_initial .= choices_current
                    end
                end
                prob = Optimization.OptimizationProblem(f, choices_initial, states, lb = choices_lb, ub = choices_ub)
                sol = solve(prob, NLopt.LN_SBPLX(), maxiters = 50, maxtime = 0.5) 
                @inbounds variables.V_4[s_i, h_i, a_i, c_i] = -sol.objective
                @inbounds variables.policy_n_4[s_i, h_i, a_i, c_i] = sol.u[1]
                @inbounds variables.policy_s_5[s_i, h_i, a_i, c_i] = sol.u[2]
            end
        end
    end
    return nothing
end