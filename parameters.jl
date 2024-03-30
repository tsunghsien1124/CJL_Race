function parameters_function(;
    #----------------------#
    # exogenous parameters #
    #----------------------#
    χ::Real=2.00,                               # CRRA coefficient
    α::Real=0.32,                               # capital share
    δ::Real=0.07,                               # depreciation rate
    σ::Real=1.0 - 1.0 / 1.44;                   # CES elasticity
    γ_0::Real=0.90;                             # share of time investment in kid's skill formation
    γ_1::Real=0.71;                             # share of time investment in kid's skill formation
    γ_2::Real=0.68;                             # share of time investment in kid's skill formation
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
    A::real=1.674,                              # TFP adjustment factor to ensure unity e_bar in equilibrium
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
    PSID_avg_earnings = 38674.2     # psidavgearn
    e_bar = 1.0                     # nmdlavgearn
    inc_bar = e_bar * (1.0 + (r / R) * (α / (1.0 - α)))

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

    # asset choices
    s_min = -g / (1.0 + r)
    s_max = e_bar * s_frac
    s_step, s_grid = h_grid_function(s_size, s_min, s_max, s_power)

    sk_min = 0.0
    sk_max = s_max * sk_frac
    sk_step, sk_grid = h_grid_function(s_size, sk_min, sk_max, s_power)

    # persistent income shock
    ϵ_MC = tauchen(ϵ_size, ρ, σ_ϵ, 0.0, 3)
    ϵ_Γ = ϵ_MC.p
    ϵ_grid = collect(ϵ_MC.state_values)
    ϵ_G = stationary_distributions(MarkovChain(ϵ_Γ, ϵ_grid))[1]

    # transitory income shock
    ν_grid, ν_Γ = adda_cooper(ν_size, 0.0, σ_ν)
    ν_Γ = ν_Γ[1, :]
    ν_G = ν_Γ

    # asset holding
    a_grid = ((range(0.0, stop=a_size - 1, length=a_size) / (a_size - 1)) .^ a_degree) * a_max

    # return values
    return (
        r=r,
        β=β,
        γ=γ,
        ρ=ρ,
        σ_ϵ=σ_ϵ,
        σ_ν=σ_ν,
        b=b,
        κ=κ,
        ψ=ψ,
        μ=μ,
        θ=θ,
        q_bar=q_bar,
        ψ_1=ψ_1,
        ψ_2=ψ_2,
        p=p,
        age_min=age_min,
        age_max=age_max,
        age_inf=age_inf,
        age_ret=age_ret,
        age_size=age_size,
        age_grid=age_grid,
        inf_grid=inf_grid,
        n_max=n_max,
        n_size=n_size,
        n_grid=n_grid,
        n_Γ=n_Γ,
        ϵ_size=ϵ_size,
        ϵ_grid=ϵ_grid,
        ϵ_Γ=ϵ_Γ,
        ϵ_G=ϵ_G,
        ν_size=ν_size,
        ν_grid=ν_grid,
        ν_Γ=ν_Γ,
        ν_G=ν_G,
        a_max=a_max,
        a_size=a_size,
        a_grid=a_grid,
        a_degree=a_degree
    )
end