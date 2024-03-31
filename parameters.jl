function parameters_function(;
    #----------------------#
    # exogenous parameters #
    #----------------------#
    χ::Real=2.00,                               # CRRA coefficient
    α::Real=0.32,                               # capital share
    δ::Real=0.07,                               # depreciation rate
    σ::Real=1.0 - 1.0 / 1.44,                   # CES elasticity
    γ_0::Real=0.90,                             # share of time investment in kid's skill formation
    γ_1::Real=0.71,                             # share of time investment in kid's skill formation
    γ_2::Real=0.68,                             # share of time investment in kid's skill formation
    r_0::Real=0.04,                             # annual interest rate
    τ_0::Real=0.10,                             # intercept of progressive income tax
    τ_1::Real=0.04,                             # slope of progressive income tax
    g::Real=0.03,                               # government welfare transfer
    q_A::Real=1.70,                             # adult equivalance scale
    p_1::Real=0.33,                             # intercept of social security payment 
    p_0::Real=0.08,                             # slope of social security payment
    τ_s::Real=0.11,                             # payroll tax
    d_0::Real=0.02,                             # public subsidy in kid's education
    d_1::Real=0.09,                             # public subsidy in kid's education
    d_2::Real=0.10,                             # public subsidy in kid's education
    κ::Real=0.19,                               # college cost
    #--------------------#
    # set in equilibrium #
    #--------------------#
    β::Real=0.98,                               # discount factor
    ν::Real=0.70,                               # weight on high school human capital
    A::Real=1.674,                              # TFP adjustment factor to ensure unity e_bar in equilibrium
    #----------------------#
    # estimated parameters #
    #----------------------#
    ϵ_μ::Vector{Float64}=[-0.13, -0.10],        # degree-dependent mean of market luck shock
    ϵ_σ::Vector{Float64}=[0.13, 0.21],          # degree-dependent standard deviation of market luck shock
    ϕ_0::Real=1.0,                              # complementarity between different stages of child investments
    ϕ_1::Real=0.0,                              # complementarity between different stages of child investments
    ϕ_2::Real=0.0,                              # complementarity between different stages of child investments
    θ::Real=0.32,                               # parental altruism
    a_ρ::Real=0.23,                             # persistence of learning ability
    a_μ::Real=0.83,                             # mean of learning ability
    a_σ::Real=0.30,                             # standard deviaiton of learning ability
    b::Real=0.83,                               # Ben-Porath human capital accumulation
    ζ::Real=3.11,                               # kid to adult human capital anchor
    ω_1::Real=0.56,                             # primary productivity of kid's investment
    ω_2::Real=0.30,                             # secondary productivity of kid's investment
    ψ_1::Real=0.23,                             # preference for kid to college with high school parent
    ψ_2::Real=0.24,                             # preference for kid to college with college parent
    #-------------------------#
    # numerical specification #
    #-------------------------#
    age_min::Integer=4,                         # min model age
    age_max::Integer=12,                        # max model age
    age_periods::Integer=6,                     # model period unit in years
    c_size::Integer=2,                          # number of academic degree 
    ϵ_size::Integer=2,                          # number of market luck shock
    a_size::Integer=3,                          # number of learning ability shock
    h_size::Integer=12,                         # number of human capital
    h_power::Integer=2,                         # grid power of human capital
    h_frac::Real=20.0,                          # upper bound scale of human capital
    hk_frac::Real=0.25,                         # upper bound scale of kid's human capital
    s_size::Integer=24,                         # number of asset choice
    s_power::Integer=2,                         # grid power of asset choice
    s_frac::Real=6.0,                           # upper bound scale of asset choices
    sk_frac::Real=1.25,                         # upper bound scale of kid's asset choices
)
    """
    contruct an immutable object of all paramters
    """
    # aggregate prices
    r = (1.0 + r_0)^age_periods - 1.0
    R = (1.0 + r_0 + δ)^age_periods - 1.0
    W_0 = (1.0 - α) * (α / R)^(α / (1.0 - α))
    W = W_0 * A

    # average earnings and income
    PSID_avg_earnings = 38674.2                         # psidavgearn
    e_bar = 1.0                                         # nmdlavgearn
    inc_bar = e_bar * (1.0 + (r / R) * (α / (1.0 - α))) # nmdlavginc

    # age
    age_size = length(age_min:age_max)

    # market luck shock
    ϵ_grid = ϵ_grid_function(ϵ_size, ϵ_μ, ϵ_σ, c_size)

    # learning ability shock
    a_Γ, a_grid, a_birth = a_grid_function(a_size, a_ρ, a_σ, a_μ)

    # human capital (parent and kid)
    h_min = 0.0
    h_max = e_bar / W * h_frac
    h_step, h_grid = h_grid_function(h_size, h_min, h_max, h_power)

    hk_max = h_max * hk_frac
    hk_step, hk_grid = h_grid_function(h_size, h_min, hk_max, h_power)

    # saving or borrowing (parent and kid)
    s_min = -g / (1.0 + r)
    s_max = e_bar * s_frac
    s_step, s_grid = h_grid_function(s_size, s_min, s_max, s_power)

    sk_min = 0.0
    sk_max = s_max * sk_frac
    sk_step, sk_grid = h_grid_function(s_size, sk_min, sk_max, s_power)

    # return values
    return (
        χ = χ,
        α = α,
        δ = δ,
        σ = σ,
        γ_0 = γ_0,
        γ_1 = γ_1,
        γ_2 = γ_2,
        r_0 = r_0,
        τ_0 = τ_0,
        τ_1 = τ_1,
        g = g,
        q_A = q_A,
        p_1 = p_1,
        p_0 = p_0,
        τ_s = τ_s,
        d_0 = d_0,
        d_1 = d_1,
        d_2 = d_2,
        κ = κ,
        β = β,
        ν = ν,
        A = A,
        ϵ_μ = ϵ_μ,
        ϵ_σ = ϵ_σ,
        ϕ_0 = ϕ_0,
        ϕ_1 = ϕ_1,
        ϕ_2 = ϕ_2,
        θ = θ,
        a_ρ = a_ρ,
        a_μ = a_μ,
        a_σ = a_σ,
        b = b,
        ζ = ζ,
        ω_1 = ω_1,
        ω_2 = ω_2,
        ψ_1 = ψ_1,
        ψ_2 = ψ_2,
        r = r,
        R = R,
        W_0 = W_0,
        W = W,
        PSID_avg_earnings = PSID_avg_earnings,
        e_bar = e_bar,
        inc_bar = inc_bar,
        age_min = age_min,
        age_max = age_max,
        age_size = age_size,
        age_periods = age_periods,
        c_size = c_size,
        ϵ_size = ϵ_size,
        ϵ_grid = ϵ_grid,
        a_size = a_size,
        a_Γ = a_Γ,
        a_grid = a_grid,
        a_birth = a_birth,
        h_size = h_size,
        h_power = h_power,
        h_frac = h_frac,
        hk_frac = hk_frac,
        h_min = h_min,
        h_max = h_max,
        h_step = h_step,
        h_grid = h_grid,
        hk_max = hk_max,
        hk_step = hk_step,
        hk_grid = hk_grid,
        s_size = s_size,
        s_power = s_power,
        s_frac = s_frac,
        sk_frac = sk_frac,
        s_min = s_min,
        s_max = s_max,
        s_step = s_step,
        s_grid = s_grid,
        sk_min = sk_min,
        sk_max = sk_max,
        sk_step = sk_step, 
        sk_grid = sk_grid,        
    )
end