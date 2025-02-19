function age_9_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    solve the financial transfer problem of households
    """
    @unpack ϵ_size, ϵ_grid, c_size, a_size, a_grid, h_size, h_grid, s_size, s_grid, s_min, s_k_min, s_k_grid = grids
    @unpack β, θ, b = parameters

    choices_candidate = collect(0.0:0.25:1.0)
    states = [0.0, 0.0, 0.0, 0.0, 0.0]
    @inbounds @views V_10_itp = linear_interpolation((s_grid, h_grid), zeros(s_size, h_size), extrapolation_bc=Interpolations.Flat())
    @inbounds @views V_4_itp = linear_interpolation(s_k_grid, zeros(s_size), extrapolation_bc=Interpolations.Flat())

    choices_initial = [0.5, 0.5, 0.5]
    choices_current = [0.5, 0.5, 0.5]
    choices_lb = [0.0, 0.0, 0.0]
    choices_ub = [1.0, 1.0, 1.0]
    function age_9_obj_function(choices::Vector{Float64}, states::Vector{Float64})
        """
        construct the objective function in age 9
            choices: (n_9, s_10, s_4)
            states: (h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i)
        """
        n_9 = choices[1]
        h_10 = states[2] * (n_9 * states[1])^b + states[1]
        earnings_9 = states[4] * states[1] * (1.0 - n_9)
        budget = max(0.0, f_function(earnings_9, states[3], 9, parameters, prices) - s_min - s_k_min)
        s_10 = choices[2] * budget
        s_4 = choices[3] * (budget - s_10)
        c_9 = budget - s_10 - s_4
        s_10 = s_10 + s_min
        s_4 = s_4 + s_k_min
        return -utility_function(c_9, parameters) - θ * V_4_itp(s_4) - β * V_10_itp(s_10, states[5] * h_10) # EV_10_itp(s_10, h_10)
    end
    f = OptimizationFunction(age_9_obj_function)

    choices_initial_h_0 = [0.5, 0.5]
    choices_current_h_0 = [0.5, 0.5]
    choices_lb_h_0 = [0.0, 0.0]
    choices_ub_h_0 = [1.0, 1.0]
    function age_9_obj_function_h_0(choices::Vector{Float64}, states::Vector{Float64})
        """
        construct the objective function in age 9
            choices: (n_9, s_10, s_4)
            states: (h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i)
        """
        n_9 = 0.0
        h_10 = states[2] * (n_9 * states[1])^b + states[1]
        earnings_9 = states[4] * states[1] * (1.0 - n_9)
        budget = max(0.0, f_function(earnings_9, states[3], 9, parameters, prices) - s_min - s_k_min)
        s_10 = choices[1] * budget
        s_4 = choices[2] * (budget - s_10)
        c_9 = budget - s_10 - s_4
        s_10 = s_10 + s_min
        s_4 = s_4 + s_k_min
        return -utility_function(c_9, parameters) - θ * V_4_itp(s_4) - β * V_10_itp(s_10, states[5] * h_10) # EV_10_itp(s_10, h_10)
    end
    f_h_0 = OptimizationFunction(age_9_obj_function_h_0)

    for c_i = 1:c_size
        @inbounds w_c = prices.w_S[c_i]
        @inbounds @views ϵ_c = ϵ_grid[:, c_i]
        ϵ_c_mean = mean(ϵ_c)
        @inbounds @views V_10_itp.itp.coefs .= variables.V_10[:, :, c_i]
        # @inbounds EV_10_itp(s_10, h_10) = (V_10_itp(s_10, ϵ_c[1] * h_10) + V_10_itp(s_10, ϵ_c[2] * h_10)) / ϵ_size
        # EV_10_itp(s_10, h_10) = V_10_itp(s_10, ϵ_c_mean * h_10)

        for c_k_i = 1:c_size, a_k_i = 1:a_size, h_k_i = 1:h_size
            @inbounds @views V_4_itp.itp.coefs .= variables.V_4[:, h_k_i, a_k_i, c_k_i]

            for a_i = 1:a_size, h_i = 1:h_size, s_i = 1:s_size
                # println("c_k_i = $c_k_i, a_k_i = $a_k_i, h_k_i = $h_k_i, a_i = $a_i, h_i = $h_i, s_i = $s_i")
                @inbounds h_9 = h_grid[h_i]
                @inbounds a_9 = a_grid[a_i]
                @inbounds s_9 = s_grid[s_i]
                states .= [h_9, a_9, s_9, w_c, ϵ_c_mean]
                if h_9 == 0.0
                    for ch_1_i in choices_candidate, ch_2_i in choices_candidate
                        choices_current_h_0 .= [ch_1_i, ch_2_i]
                        if age_9_obj_function_h_0(choices_current_h_0, states) < age_9_obj_function_h_0(choices_initial_h_0, states)
                            choices_initial_h_0 .= choices_current_h_0
                        end
                    end
                    prob = Optimization.OptimizationProblem(f_h_0, choices_initial_h_0, states, lb = choices_lb_h_0, ub = choices_ub_h_0)
                    sol = solve(prob, NLopt.LN_SBPLX(), maxiters = 50, maxtime = 0.5)
                    @inbounds variables.V_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = -sol.objective
                    @inbounds variables.policy_n_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = 0.0
                    @inbounds variables.policy_s_10[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = sol.u[1]
                    @inbounds variables.policy_s_4[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = sol.u[2]
                    # for ch_1_i in choices_candidate, ch_2_i in choices_candidate, ch_3_i in choices_candidate
                    #     choices_current .= [ch_1_i, ch_2_i, ch_3_i]
                    #     if age_9_obj_function(choices_current, states) < age_9_obj_function(choices_initial, states)
                    #     choices_initial .= choices_current
                    #     end
                    # end
                    # choices_initial[1] = 0.0 
                    # prob = Optimization.OptimizationProblem(f, choices_initial, states, lb = choices_lb, ub = [0.0, 1.0, 1.0])
                    # sol = solve(prob, NLopt.LN_BOBYQA(), maxiters = 50, maxtime = 0.5)
                    # @inbounds variables.V_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = -sol.objective
                    # @inbounds variables.policy_n_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = sol.u[1]
                    # @inbounds variables.policy_s_10[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = sol.u[2]
                    # @inbounds variables.policy_s_4[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = sol.u[3]
                else
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
    end
    return nothing
end

# function age_9_function!(V_10::Array{Float64,3}, V_9::Array{Float64,7}, V_4::Array{Float64,4}, policy_n_9::Array{Float64,7}, policy_s_10::Array{Float64,7}, policy_s_4::Array{Float64,7}, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
#     """
#     solve the financial transfer problem of households
#     """
#     @unpack ϵ_size, ϵ_grid, c_size, a_size, a_grid, h_size, h_grid, s_size, s_grid, s_min, s_k_min, s_k_grid = grids
#     @unpack β, θ, b = parameters

#     choices_initial = [0.5, 0.5, 0.5]
#     choices_lb = [0.0, 0.0, 0.0]
#     choices_ub = [1.0, 1.0, 1.0]
#     @inbounds @views V_10_itp = linear_interpolation((s_grid, h_grid), zeros(s_size, h_size), extrapolation_bc=Line())
#     @inbounds @views V_4_itp = linear_interpolation(s_k_grid, zeros(s_size), extrapolation_bc=Line())

#     for c_i = 1:c_size
#         @inbounds w_c = prices.w_S[c_i]
#         @inbounds @views ϵ_c = ϵ_grid[:, c_i]
#         @inbounds @views V_10_itp.itp.coefs .= V_10[:, :, c_i]
#         @inbounds EV_10_itp(s_10, h_10) = (V_10_itp(s_10, ϵ_c[1] * h_10) + V_10_itp(s_10, ϵ_c[2] * h_10)) / ϵ_size
#         for c_k_i = 1:c_size, a_k_i = 1:a_size, h_k_i = 1:h_size
#             @inbounds @views V_4_itp.itp.coefs .= V_4[:, h_k_i, a_k_i, c_k_i]
#             for a_i = 1:a_size, h_i = 1:h_size, s_i = 1:s_size
#                 @inbounds h_9 = h_grid[h_i]
#                 @inbounds a_9 = a_grid[a_i]
#                 @inbounds s_9 = s_grid[s_i]
#                 function age_9_obj_function(choices::Vector{Float64})
#                     """
#                     construct the objective function in age 9
#                         choices: (n_9, s_10, s_4)
#                         states: (h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i)
#                     """
#                     n_9 = choices[1]
#                     h_10 = a_9 * (n_9 * h_9)^b + h_9
#                     earnings_9 = w_c * h_9 * (1.0 - n_9)
#                     budget = max(0.0, f_function(earnings_9, s_9, 9, parameters, prices) - s_min - s_k_min)
#                     s_10 = choices[2] * budget
#                     s_4 = choices[3] * (budget - s_10)
#                     c_9 = budget - s_10 - s_4
#                     s_10 = s_10 + s_min
#                     s_4 = s_4 + s_k_min
#                     return -utility_function(c_9, parameters) - θ * V_4_itp(s_4) - β * EV_10_itp(s_10, h_10)
#                 end
#                 res = optimize(age_9_obj_function, choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()))
#                 @inbounds V_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = -res.minimum
#                 @inbounds policy_n_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[1]
#                 @inbounds policy_s_10[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[2]
#                 @inbounds policy_s_4[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[3]
#             end
#         end
#     end
#     return nothing
# end