function age_9_function_th!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple) # V_10_temp::Array{Float64, 9},
    """
    solve the financial transfer problem of households
    """
    @unpack ϵ_size, ϵ_grid, c_size, a_size, a_grid, h_size, h_grid, s_size, s_grid, s_min, s_k_min, loop_age_10 = grids
    @unpack β, θ, b = parameters
    choices_initial = [0.5, 0.5, 0.5]
    choices_lb = [0.0, 0.0, 0.0]
    choices_ub = [1.0, 1.0, 1.0]
    # V_10_temp = zeros(ϵ_size, 2)
    # h_9, a_9, s_9, w_c, h_10 = 0.0, 0.0, 0.0, 0.0, 0.0

    # V_9 = zeros(h_size, a_size, c_size, s_size, h_size, a_size, c_size)
    # Threads.@threads 
    # @batch 
    for (c_i, a_i, h_i, s_i, c_k_i, a_k_i, h_k_i) in loop_age_10

        h_9 = h_grid[h_i]
        a_9 = a_grid[a_i]
        s_9 = s_grid[s_i]
        w_c = prices.w_S[c_i]
        # V_10_temp = zeros(ϵ_size, 2)

        function age_9_obj_function(choices::Vector{Float64})
            """
            construct the objective function in age 9
                choices: (n_9, s_10, s_4)
                states: (h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i)
            """
            # @views @inbounds V_10_temp[:, :, h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] .= 0.0
        
            # results container
            # V_10_temp = zeros(ϵ_size, 2)
            V_10_temp = 0.0
            s_k_ind_lower, s_k_wgt_lower = 0, 0.0
            s_ind_lower, s_wgt_lower = 0, 0.0
        
            # compute implied value for the given choices
            n_9 = choices[1]
            h_10 = a_9 * (n_9 * h_9)^b + h_9
            # earnings_9 = w_c * h_9 * (1.0 - n_9)
            budget = max(0.0, f_function(w_c * h_9 * (1.0 - n_9), s_9, 9, parameters, prices) - s_min - s_k_min)
            s_10 = choices[2] * budget
            s_4 = choices[3] * (budget - s_10)
            c_9 = budget - s_10 - s_4
            s_10 = s_10 + s_min
            s_4 = s_4 + s_k_min
        
            # interpolate V_4
            s_k_ind_lower, s_k_wgt_lower = locate_s_kid_function(s_4, grids)
            # locate_s_kid_function!(s_k_ind_lower, s_k_wgt_lower, s_4, grids)
            @inbounds itp_V_4 = s_k_wgt_lower * variables.V_4[s_k_ind_lower, h_k_i, a_k_i, c_k_i] + (1.0 - s_k_wgt_lower) * variables.V_4[s_k_ind_lower+1, h_k_i, a_k_i, c_k_i]
        
            # interpolate V_10
            # s_ind_lower, s_wgt_lower = locate_s_function(s_10, grids)
            locate_s_function!(s_ind_lower, s_wgt_lower, s_10, grids)
            for ϵ_i = 1:ϵ_size
                # h_10_ϵ = ϵ_grid[ϵ_i, c_i] * h_10
                @inbounds h_ind_lower, h_wgt_lower = locate_h_function(ϵ_grid[ϵ_i, c_i] * h_10, grids)
                # @inbounds V_10_temp[ϵ_i, 1] = s_wgt_lower * variables.V_10[s_ind_lower, h_ind_lower, c_i] + (1.0 - s_wgt_lower) * variables.V_10[s_ind_lower+1, h_ind_lower, c_i]
                # @inbounds V_10_temp[ϵ_i, 2] = s_wgt_lower * variables.V_10[s_ind_lower, h_ind_lower, c_i] + (1.0 - s_wgt_lower) * variables.V_10[s_ind_lower+1, h_ind_lower+1, c_i]
                # @inbounds V_10_temp[ϵ_i, 1] = h_wgt_lower * V_10_temp[ϵ_i, 1] + (1.0 - h_wgt_lower) * V_10_temp[ϵ_i, 2]
                # @views @inbounds V_10_temp[ϵ_i, 1, h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = s_wgt_lower * variables.V_10[s_ind_lower, h_ind_lower, c_i] + (1.0 - s_wgt_lower) * variables.V_10[s_ind_lower+1, h_ind_lower, c_i]
                # @views @inbounds V_10_temp[ϵ_i, 2, h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = s_wgt_lower * variables.V_10[s_ind_lower, h_ind_lower, c_i] + (1.0 - s_wgt_lower) * variables.V_10[s_ind_lower+1, h_ind_lower+1, c_i]
                # @views @inbounds V_10_temp[ϵ_i, 1, h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = h_wgt_lower * V_10_temp[ϵ_i, 1] + (1.0 - h_wgt_lower) * V_10_temp[ϵ_i, 2]
                @inbounds V_10_temp += h_wgt_lower * (s_wgt_lower * variables.V_10[s_ind_lower, h_ind_lower, c_i] + (1.0 - s_wgt_lower) * variables.V_10[s_ind_lower+1, h_ind_lower, c_i])
                @inbounds V_10_temp += (1.0 - h_wgt_lower) * (s_wgt_lower * variables.V_10[s_ind_lower, h_ind_lower, c_i] + (1.0 - s_wgt_lower) * variables.V_10[s_ind_lower+1, h_ind_lower+1, c_i])
            end
            # @views @inbounds itp_V_10 = sum(V_10_temp[:, 1, h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i]) / ϵ_size
            @views @inbounds itp_V_10 = V_10_temp / ϵ_size
        
            # return implied V_9 for the given choices
            return -utility_function(c_9, parameters) - θ * itp_V_4 - β * itp_V_10
        end

        res = optimize(choices -> age_9_obj_function(choices), choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()))
        @inbounds variables.V_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = -res.minimum
        @inbounds variables.policy_n_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[1]
        @inbounds variables.policy_s_10[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[2]
        @inbounds variables.policy_s_4[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[3]
        # res_obj = OptimizationProblem(age_9_obj_function, choices_initial, lb = choices_lb, ub = choices_ub)
        # res_sol = solve(res_obj, BBO_adaptive_de_rand_1_bin_radiuslimited())
        # @inbounds variables.V_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = -res_sol.objective
        # @inbounds variables.policy_n_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res_sol.u[1]
        # @inbounds variables.policy_s_10[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res_sol.u[2]
        # @inbounds variables.policy_s_4[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res_sol.u[3]
        # println(states)
    end
end
