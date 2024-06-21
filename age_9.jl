# function age_9_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
#     """
#     solve the financial transfer problem of households
#     """
#     @unpack c_size, a_size, h_size, h_grid, s_size = grids
#     choices_initial = [0.5, 0.5, 0.5]
#     choices_lb = [0.0, 0.0, 0.0]
#     choices_ub = [1.0, 1.0, 1.0]
    
#     # V_9 = zeros(h_size, a_size, c_size, s_size, h_size, a_size, c_size)
#     for c_i = 1:c_size, a_i = 1:a_size, h_i = 1:h_size, s_i = 1:s_size, c_k_i = 1:c_size, a_k_i = 1:a_size, h_k_i = 1:h_size
#         states = [h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i]
#         # choice_initial_step = 0.25
#         # choices_candidate = [(x_1, x_2, x_3) for x_1 = 0.0:choice_initial_step:1.0, x_2 = 0.0:choice_initial_step:1.0, x_3 = 0.0:choice_initial_step:1.0]
#         # choices_initial = collect(choices_candidate[findmin([age_9_obj_function([x_1, x_2, x_3], states, variables, prices, parameters) for x_1 = 0.0:choice_initial_step:1.0, x_2 = 0.0:choice_initial_step:1.0, x_3 = 0.0:choice_initial_step:1.0])[2]])
#         # # choices_initial = collect(choices_candidate[findmin(age_9_obj_function.(choices_candidate))[2]])
#         # choices_initial[choices_initial .== 0.0] .= 0.1
#         # choices_initial[choices_initial .== 1.0] .= 0.9
#         res = optimize(choices -> age_9_obj_function(choices, states, variables, prices, parameters, grids), choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()));
#         @inbounds variables.V_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = -res.minimum
#         @inbounds variables.policy_n_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[1]
#         @inbounds variables.policy_s_10[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[2]
#         @inbounds variables.policy_s_4[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[3]
#         # println(states)
#     end
# end

function age_9_function_th!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    solve the financial transfer problem of households
    """
    @timeit tmr "unpack" @unpack c_size, a_size, h_size, h_grid, s_size = grids
    @timeit tmr "int" begin
        choices_initial = [0.5, 0.5, 0.5]
        choices_lb = [0.0, 0.0, 0.0]
        choices_ub = [1.0, 1.0, 1.0]
    end

    # V_9 = zeros(h_size, a_size, c_size, s_size, h_size, a_size, c_size)
    @timeit tmr "loop" begin
        # Threads.@threads
        for (c_i, a_i, h_i, s_i, c_k_i, a_k_i, h_k_i) in collect(Iterators.product(1:c_size, 1:a_size, 1:h_size, 1:s_size, 1:c_size, 1:a_size, 1:h_size))
            states = [h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i]
            # choice_initial_step = 0.25
            # choices_candidate = [(x_1, x_2, x_3) for x_1 = 0.0:choice_initial_step:1.0, x_2 = 0.0:choice_initial_step:1.0, x_3 = 0.0:choice_initial_step:1.0]
            # choices_initial = collect(choices_candidate[findmin([age_9_obj_function([x_1, x_2, x_3], states, variables, prices, parameters) for x_1 = 0.0:choice_initial_step:1.0, x_2 = 0.0:choice_initial_step:1.0, x_3 = 0.0:choice_initial_step:1.0])[2]])
            # # choices_initial = collect(choices_candidate[findmin(age_9_obj_function.(choices_candidate))[2]])
            # choices_initial[choices_initial .== 0.0] .= 0.1
            # choices_initial[choices_initial .== 1.0] .= 0.9
            @timeit tmr "opt" begin
                res = optimize(choices -> age_9_obj_function(choices, states, variables, prices, parameters, grids), choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()))
                # res = optimize(choices -> age_9_obj_function(choices), choices_lb, choices_ub, choices_initial, Fminbox(NelderMead()))
            end
            @inbounds variables.V_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = -res.minimum
            @inbounds variables.policy_n_9[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[1]
            @inbounds variables.policy_s_10[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[2]
            @inbounds variables.policy_s_4[h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i] = res.minimizer[3]
            # println(states)
        end
    end
end

# function age_9_obj_function(choices)
#     """
#     construct the objective function in age 9
#         choices: (n_9, s_10, s_4)
#         states: (h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i)
#     """
#     @unpack h_grid, a_grid, s_grid, β, θ, b, s_min, s_k_min = parameters

#     # check if the number of choices is correct
#     if length(choices) != 3
#         error("dimension mismatch!")
#     end

#     # compute implied value for the given choices
#     h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i = states
#     n_9 = choices[1]
#     h_9 = h_grid[h_i]
#     a_9 = a_grid[a_i]
#     s_9 = s_grid[s_i]
#     earnings_9 = prices.w_S[c_i] * h_9 * (1.0 - n_9)
#     h_10 = a_9 * (n_9 * h_9)^b + h_9
#     budget = max(0.0, f_function(earnings_9, s_9, 9, parameters) - s_min - s_k_min)
#     s_10 = choices[2] * budget
#     s_4 = choices[3] * (budget - s_10)
#     c_9 = budget - s_10 - s_4
#     s_10 = s_10 + s_min
#     s_4 = s_4 + s_k_min
#     s_k_ind, s_k_wgt = locate_s_kid_function(s_4, parameters)
#     itp_V_4 = sum(s_k_wgt .* variables.V_4[s_k_ind[1]:s_k_ind[2], h_k_i, a_k_i, c_k_i])
#     itp_V_10 = V_10_itp_function(s_10, h_10, c_i, variables, parameters)
#     return -utility_function(c_9, parameters) - θ * itp_V_4 - β * itp_V_10
# end

@timeit tmr "obj" function age_9_obj_function(choices::Vector{Float64}, states::Vector{Int64}, variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple)
    """
    construct the objective function in age 9
        choices: (n_9, s_10, s_4)
        states: (h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i)
    """
    @unpack β, θ, b = parameters
    @unpack h_grid, a_grid, s_grid, s_min, s_k_min = grids

    # check if the number of choices is correct
    if length(choices) != 3
        error("dimension mismatch!")
    end

    # compute implied value for the given choices
    h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i = states
    n_9 = choices[1]
    h_9 = h_grid[h_i]
    a_9 = a_grid[a_i]
    s_9 = s_grid[s_i]
    earnings_9 = prices.w_S[c_i] * h_9 * (1.0 - n_9)
    h_10 = a_9 * (n_9 * h_9)^b + h_9
    @timeit tmr "budget" budget = max(0.0, f_function(earnings_9, s_9, 9, parameters, prices) - s_min - s_k_min)
    s_10 = choices[2] * budget
    s_4 = choices[3] * (budget - s_10)
    c_9 = budget - s_10 - s_4
    s_10 = s_10 + s_min
    s_4 = s_4 + s_k_min
    @timeit tmr "loc_s_kid" s_k_ind, s_k_wgt = locate_s_kid_function(s_4, grids)
    itp_V_4 = sum(s_k_wgt .* variables.V_4[s_k_ind[1]:s_k_ind[2], h_k_i, a_k_i, c_k_i])
    itp_V_10 = V_10_itp_function(s_10, h_10, c_i, variables, grids)
    return -utility_function(c_9, parameters) - θ * itp_V_4 - β * itp_V_10
end

function age_9_obj_function(choices::Vector{Float64})
    """
    construct the objective function in age 9
        choices: (n_9, s_10, s_4)
        states: (h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i)
    """
    @unpack β, θ, b = parameters
    @unpack h_grid, a_grid, s_grid, s_min, s_k_min = grids

    # check if the number of choices is correct
    if length(choices) != 3
        error("dimension mismatch!")
    end

    # compute implied value for the given choices
    h_k_i, a_k_i, c_k_i, s_i, h_i, a_i, c_i = states
    n_9 = choices[1]
    h_9 = h_grid[h_i]
    a_9 = a_grid[a_i]
    s_9 = s_grid[s_i]
    earnings_9 = prices.w_S[c_i] * h_9 * (1.0 - n_9)
    h_10 = a_9 * (n_9 * h_9)^b + h_9
    budget = max(0.0, f_function(earnings_9, s_9, 9, parameters, prices) - s_min - s_k_min)
    s_10 = choices[2] * budget
    s_4 = choices[3] * (budget - s_10)
    c_9 = budget - s_10 - s_4
    s_10 = s_10 + s_min
    s_4 = s_4 + s_k_min
    s_k_ind, s_k_wgt = locate_s_kid_function(s_4, grids)
    itp_V_4 = sum(s_k_wgt .* variables.V_4[s_k_ind[1]:s_k_ind[2], h_k_i, a_k_i, c_k_i])
    itp_V_10 = V_10_itp_function(s_10, h_10, c_i, variables, grids)
    return -utility_function(c_9, parameters) - θ * itp_V_4 - β * itp_V_10
end


@timeit tmr "v_10" function V_10_itp_function(s_10::Real, h_10::Real, c_i::Integer, variables::Mutable_Variables, grids::NamedTuple)
    """
    interpolate value function at age 10
        s_10: savings made at age 9
        h_10: imputed human capital at age 9 prior to market luck shock
        c_i: education indicator
    """
    @unpack ϵ_size, ϵ_grid = grids
    @timeit tmr "V_10_temp" V_10_temp = zeros(ϵ_size, 2)
    @timeit tmr "loc_s_10" s_ind, s_wgt = locate_s_function(s_10, grids)
    @timeit tmr "loc_h_10_ϵ" begin
    for ϵ_i = 1:ϵ_size
        h_10_ϵ = ϵ_grid[ϵ_i, c_i] * h_10
        h_ind, h_wgt = locate_h_function(h_10_ϵ, grids)
        V_10_temp[ϵ_i, 1] = sum(s_wgt .* variables.V_10[s_ind[1]:s_ind[2], h_ind[1], c_i])
        V_10_temp[ϵ_i, 2] = sum(s_wgt .* variables.V_10[s_ind[1]:s_ind[2], h_ind[2], c_i])
        V_10_temp[ϵ_i, 1] = sum(h_wgt .* V_10_temp[ϵ_i, :])
    end
    end
    return sum(V_10_temp[:, 1]) / ϵ_size
end