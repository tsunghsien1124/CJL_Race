function ϵ_grid_function(ϵ_size::Integer, ϵ_μ::Vector{Float64}, ϵ_σ::Vector{Float64}, c_size::Integer)
    """
    generate ϵ grid using Kennan (2006)'s Proposition B
        ϵ_size: grid size of market luck shock
        ϵ_μ: degree-dependent mean of market luck shock
        ϵ_σ: degree-dependent standard deviation of market luck shock
        c_size: grid size of degree indicator (1: high school; 2: college), i.e., c_size = 2 
    """
    ϵ_grid = zeros(ϵ_size, c_size)
    for c_i = 1:c_size
        dist = Normal(ϵ_μ[c_i], ϵ_σ[c_i])

        for ϵ_i = 1:ϵ_size
            ϵ_grid[ϵ_i, c_i] = exp(quantile(dist, (ϵ_i - 0.5) / ϵ_size))
        end
    end
    return ϵ_grid
end

function a_grid_function(a_size::Integer, a_ρ::Real, a_σ::Real, a_μ::Real)
    """
    generate a grid using Rouwenhorst from QuantEcon
        a_size: grid size of ability shock
        a_ρ: persistence of ability shock
        a_σ: standard deviaiton of ability shock
        a_μ: mean of ability shock
    """
    # (1) using QuantEcon's function "rouwenhorst"
    # a_MC = rouwenhorst(a_size, a_ρ, a_σ, a_μ)
    # a_Γ = a_MC.p
    # a_grid = collect(a_MC.state_values)
    # a_birth = stationary_distributions(MarkovChain(a_Γ, a_grid))[1]
    # return a_Γ, a_grid, a_birth

    # (2) following Lee and Seshadri's code implementation
    a_MC = rouwenhorst(a_size, a_ρ, a_σ)
    a_Γ = a_MC.p
    a_grid = exp.(collect(a_MC.state_values))
    a_birth = stationary_distributions(MarkovChain(a_Γ, a_grid))[1]
    a_grid = a_grid ./ sum(a_birth .* a_grid)
    a_grid = a_grid .* a_μ
    # mean_log_a = sum(a_birth .* log.(a_grid)) # mean of log normalized a grid
    return a_Γ, a_grid, a_birth
end

function h_grid_function(h_size::Integer, h_min::Real, h_max::Real, h_power::Real)
    """
    generate h power-grid
        h_size: grid size of human capital
        h_min: minimum of human capital
        h_max: maximum of human capital
        h_power: power of human capital grid
    """
    h_step = (h_max - h_min) / (h_size - 1)^h_power
    h_grid = zeros(h_size)
    for h_i = 1:h_size
        h_grid[h_i] = h_step * (h_i - 1)^h_power + h_min
    end
    return h_step, h_grid
end

function h_k_grid_function(h_k_size::Integer, h_k_min::Real, h_k_max::Real, h_k_power::Real)
    """
    generate h_k power-grid
        h_k_size: grid size of kid's human capital
        h_k_min: minimum of kid's human capital
        h_k_max: maximum of kid's human capital
        h_k_power: power of kid's human capital grid
    """
    h_k_step = (h_k_max - h_k_min) / (h_k_size - 1)^h_k_power
    h_k_grid = zeros(h_k_size)
    for h_k_i = 1:h_k_size
        h_k_grid[h_k_i] = h_k_step * (h_k_i - 1)^h_k_power + h_k_min
    end
    return h_k_step, h_k_grid
end

function s_grid_function(s_size::Integer, s_min::Real, s_max::Real, s_power::Real)
    """
    generate s power-grid
        s_size: grid size of savings
        s_min: minimum of savings
        s_max: maximum of savings
        s_power: power of savings grid
    """
    s_step = (s_max - s_min) / (s_size - 1)^s_power
    s_grid = zeros(s_size)
    for s_i = 1:s_size
        s_grid[s_i] = s_step * (s_i - 1)^s_power + s_min
    end
    return s_step, s_grid
end

function s_k_grid_function(s_k_size::Integer, s_k_min::Real, s_k_max::Real, s_k_power::Real)
    """
    generate s power-grid
        s_k_size: grid size of kid's savings
        s_k_min: minimum of kid's savings
        s_k_max: maximum of kid's savings
        s_k_power: power of kid's savings grid
    """
    s_k_step = (s_k_max - s_k_min) / (s_k_size - 1)^s_k_power
    s_k_grid = zeros(s_k_size)
    for s_k_i = 1:s_k_size
        s_k_grid[s_i] = s_k_step * (s_k_i - 1)^s_k_power + s_k_min
    end
    return s_k_step, s_k_grid
end

function print_grid_function(grid::Vector{Float64}, filename::String)
    """
    print out any grid and save it as txt file with associated grid name
    """
    open(pwd() * "\\" * filename * ".txt", "w") do io
        writedlm(io, grid)
    end
end

function locate_h_function(h::Real, grids::NamedTuple)
    """
    locate h on its grid with indices and convex weights
        h: human capital
        parameters: collection of parameters
    """
    @unpack h_size, h_step, h_min, h_power, h_grid = grids
    h_ind_lower = max(1, min(h_size - 1, floor(Int, ((h - h_min) / h_step)^(1.0 / h_power) + 1)))
    h_ind_upper = h_ind_lower + 1
    h_wgt_lower = max(0.0, min(1.0, (h_grid[h_ind_upper] - h) / (h_grid[h_ind_upper] - h_grid[h_ind_lower])))
    # h_wgt_upper = 1.0 - h_wgt_lower
    return h_ind_lower, h_wgt_lower #, h_ind_upper, h_wgt_upper
end

function locate_h_kid_function(h_k::Real, grids::NamedTuple)
    """
    locate h_k on its grid with indices and convex weights
        h_k: kid's human capital
        parameters: collection of parameters
    """
    @unpack h_size, h_k_step, h_min, h_power, h_k_grid = grids
    h_k_ind_lower = max(1, min(h_size - 1, floor(Int, ((h_k - h_min) / h_k_step)^(1.0 / h_power) + 1)))
    h_k_ind_upper = h_k_ind_lower + 1
    h_k_wgt_lower = (h_k_grid[h_k_ind_upper] - h_k) / (h_k_grid[h_k_ind_upper] - h_k_grid[h_k_ind_lower])
    h_k_wgt_upper = 1.0 - h_k_wgt_lower
    h_k_ind = [h_k_ind_lower, h_k_ind_upper]
    h_k_wgt = [h_k_wgt_lower, h_k_wgt_upper]
    h_k_wgt = max.(0.0, min.(1.0, h_k_wgt))
    return h_k_ind, h_k_wgt
end

function locate_s_function(s::Real, grids::NamedTuple)
    """
    locate s on its grid with indices and convex weights
        s: savings
        parameters: collection of parameters
    """
    @unpack s_size, s_step, s_min, s_power, s_grid = grids
    s_ind_lower = max(1, min(s_size - 1, floor(Int, ((s - s_min) / s_step)^(1.0 / s_power) + 1)))
    s_ind_upper = s_ind_lower + 1
    s_wgt_lower = max(0.0, min(1.0, (s_grid[s_ind_upper] - s) / (s_grid[s_ind_upper] - s_grid[s_ind_lower])))
    # s_wgt_upper = 1.0 - s_wgt_lower
    return s_ind_lower, s_wgt_lower #, s_ind_upper, s_wgt_upper
end

function locate_s_kid_function(s_k, grids::NamedTuple)
    """
    locate s_k on its grid with indices and convex weights
        s_k: kid's savings
        parameters: collection of parameters
    """
    @unpack s_size, s_k_step, s_k_min, s_power, s_k_grid = grids
    s_k_ind_lower = max(1, min(s_size - 1, floor(Int, ((s_k - s_k_min) / s_k_step)^(1.0 / s_power) + 1))) # findlast(s_k .>= s_k_grid)
    s_k_ind_upper = s_k_ind_lower + 1
    s_k_wgt_lower = max(0.0, min(1.0, (s_k_grid[s_k_ind_upper] - s_k) / (s_k_grid[s_k_ind_upper] - s_k_grid[s_k_ind_lower])))
    # s_k_wgt_upper = 1.0 - s_k_wgt_lower
    return s_k_ind_lower, s_k_wgt_lower #, s_k_ind_upper, s_k_wgt_upper
end

function utility_function(c::Real, parameters::NamedTuple)
    """
    utility function with normalization factor
        c: comsumption
        χ: CRRA coefficient
    """
    @unpack β, θ, χ = parameters
    if c > 0.0
        normal_factor = (1.0 - β) * (1 - β^5 * θ) / (1.0 - β^9)
        return normal_factor * (c^(1.0 - χ) / (1.0 - χ))
    else
        return -10.0^9
    end
end

function tax_rate_function(y::Real, parameters::NamedTuple, prices::Mutable_Prices)
    """
    tax rate relative to mean income y_bar = e_bar + r*s_bar
        y: period income
        parameters: collection of parameters
    """
    @unpack τ_0, τ_1 = parameters
    if y <= 0.0
        return 0.0
    else
        tax_rate = τ_0 + τ_1 * log(y / prices.inc_bar)
        return max(0.0, min(1.0, tax_rate))
    end
end

function f_function(e::Real, s::Real, j::Integer, parameters::NamedTuple, prices::Mutable_Prices)
    """
    after-tax income
        e: earnings
        s: savings
        j: age
        parameters: collection of parameters
    """
    @unpack τ_s, q_A, g = parameters
    y = e + prices.r * s
    after_tax_income = (1.0 - tax_rate_function(y, parameters, prices)) * y - (1.0 - τ_s) * e
    if (j == 5) || (j == 6) || (j == 7)
        return max(0.0, after_tax_income) + q_A * g
    else
        return max(0.0, after_tax_income) + g
    end
end
