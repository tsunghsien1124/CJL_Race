function age_10_function!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple)
    """
    solve the retirment problem of households
    """
    @unpack c_size, h_size, h_grid, s_size, s_grid, β, τ_k, r, p_0, p_1, g, χ = parameters

    r_tilde = (1 - τ_k) * r
    aux_coeff = 1.0 + β^(1.0 / χ) * (1.0 + r_tilde)^(1.0 / (χ - 1.0))
    aux_coeff = 1.0 + β^(1.0 / χ) * (1.0 + r_tilde)^(1.0 / (χ - 1.0)) * aux_coeff
    for c_i = 1:c_size, h_i = 1:h_size
        earnings = prices.w_S[c_i] * h_grid[h_i]
        for s_i = 1:s_size
            total_income = (2.0 + r_tilde / (1.0 + r_tilde)^2.0) * (p_0 + p_1 * earnings + g)
            total_income = f_function(earnings, s_grid[s_i], 10, parameters) + s_grid[s_i] + total_income
            c_10 = total_income / aux_coeff
            variables.V_10[s_i, h_i, c_i] = utility_function(c_10, parameters) * aux_coeff
        end
    end
end