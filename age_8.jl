function age_8_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    solve the household problem for kid's college decision
    """
    @unpack ψ_1, ψ_2 = parameters
    @unpack c_size, a_size, h_size, h_grid, s_size = grids

    age_8_W_function!(variables, prices, parameters, grids)

    for c_i = 1:c_size, a_i = 1:a_size, h_i = 1:h_size, s_i = 1:s_size, a_k_i = 1:a_size, h_k_i = 1:h_size
        ψ_aux = c_i == 1 ? ψ_1 : ψ_2
        if variables.W_8[h_k_i, a_k_i, 1, s_i, h_i, a_i, c_i] < variables.W_8[h_k_i, a_k_i, 2, s_i, h_i, a_i, c_i] + ψ_aux
            @inbounds variables.V_8[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = variables.W_8[h_k_i, a_k_i, 2, s_i, h_i, a_i, c_i] + ψ_aux
            @inbounds variables.policy_c_3[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = 1.0
        else
            @inbounds variables.V_8[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = variables.W_8[h_k_i, a_k_i, 1, s_i, h_i, a_i, c_i]
            @inbounds variables.policy_c_3[h_k_i, a_k_i, s_i, h_i, a_i, c_i] = 0.0
        end
    end
end

function age_8_W_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    solve the household problem given kid's college state
    """
    @unpack c_size, a_size, a_grid, h_size, h_grid, h_k_grid, s_size, s_grid, s_min, ϵ_size, ϵ_grid = grids
    @unpack b, q, β, κ_tilde = parameters

    choices_candidate = collect(0.1:0.2:0.9)
    states = zeros(9) # h_8, a_8, s_8, w_c, ϵ_c_mean, h_k_3, a_k_3, w_c_k, ϵ_c_k_mean
    @inbounds @views V_9_itp = linear_interpolation((h_k_grid, s_grid, h_grid), zeros(h_size, s_size, h_size), extrapolation_bc=Interpolations.Flat())

    choices_initial = [0.5, 0.5, 0.5] # n_8, s_9, n_k_3
    choices_current = [0.5, 0.5, 0.5]
    choices_lb = [0.0, 0.0, 0.0]
    choices_ub = [1.0, 1.0, 1.0]
    function age_8_obj_function(choices::Vector{Float64}, states::Vector{Float64})
        """
        construct the objective function in age 8
            choices: (n_8, s_9, n_k_3)
            states: (h_k_3, a_k, c_k, s_8, h_8, a, c)
        """
        h_8, a_8, s_8, w_c, ϵ_c_mean, h_k_3, a_k_3, w_c_k, ϵ_c_k_mean = states
        n_8 = choices[1]
        h_9 = a_8 * (n_8 * h_8)^b + h_8
        earnings_8 = w_c * h_8 * (1.0 - n_8)
        n_k_3 = w_c_k == prices.w_S[2] ? κ_tilde + (1.0 - κ_tilde) * choices[3] : choices[3]
        h_k_4 = a_k_3 * (n_k_3 * h_k_3)^b + h_k_3
        earnings_k_3 = w_c_k * h_k_3 * (1.0 - n_k_3)
        budget = max(0.0, f_function(earnings_8, s_8, 8, parameters, prices) + f_function(earnings_k_3, 0.0, 3, parameters, prices) - s_min)
        s_9 = choices[2] * budget
        c_8 = budget - s_9
        s_9 = s_9 + s_min
        return -q * utility_function(c_8, parameters) - β * V_9_itp(ϵ_c_k_mean * h_k_4, s_9, ϵ_c_mean * h_9)
    end
    f = OptimizationFunction(age_8_obj_function)

    for c_i = 1:c_size, c_k_i = 1:c_size
        @inbounds w_c = prices.w_S[c_i]
        @inbounds @views ϵ_c = ϵ_grid[:, c_i]
        ϵ_c_mean = mean(ϵ_c)
        @inbounds w_c_k = prices.w_S[c_k_i]
        @inbounds @views ϵ_c_k = ϵ_grid[:, c_k_i]
        ϵ_c_k_mean = mean(ϵ_c_k)
        for a_i = 1:a_size, a_k_i = 1:a_size
            @inbounds a_8 = a_grid[a_i]
            @inbounds a_k_3 = a_grid[a_k_i]
            @inbounds @views V_9_itp.itp.coefs .= variables.V_9[:, a_k_i, c_k_i, :, :, a_i, c_i]
            for h_i = 1:h_size, s_i = 1:s_size, h_k_i = 1:h_size
                @inbounds h_8 = h_grid[h_i]
                @inbounds s_8 = s_grid[s_i]
                @inbounds h_k_3 = h_k_grid[h_k_i]
                states .= [h_8, a_8, s_8, w_c, ϵ_c_mean, h_k_3, a_k_3, w_c_k, ϵ_c_k_mean]
                for ch_1_i in choices_candidate, ch_2_i in choices_candidate, ch_3_i in choices_candidate
                    choices_current .= [ch_1_i, ch_2_i, ch_3_i]
                    if age_8_obj_function(choices_current, states) < age_8_obj_function(choices_initial, states)
                    choices_initial .= choices_current
                    end
                end
                prob = Optimization.OptimizationProblem(f, choices_initial, states, lb = choices_lb, ub = choices_ub)
                sol = solve(prob, NLopt.LN_SBPLX(), maxiters = 100, maxtime = 0.5)
                @inbounds variables.W_8[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = -sol.objective
                @inbounds variables.policy_n_8[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = sol.u[1]
                @inbounds variables.policy_s_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = sol.u[2]
                @inbounds variables.policy_n_3[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = sol.u[3]
            end
        end
    end
end
