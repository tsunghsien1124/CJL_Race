function parameters_function(;
    #----------------------#
    # exogenous parameters #
    #----------------------#
    age_periods::Integer=6,                     # model period unit in years
    χ::Real=2.00,                               # CRRA coefficient
    α::Real=0.32,                               # capital share
    δ::Real=0.07,                               # depreciation rate
    σ::Real=1.0 - 1.0 / 1.44,                   # CES elasticity
    γ_0::Real=0.90,                             # share of time investment in kid's s_kill formation
    γ_1::Real=0.71,                             # share of time investment in kid's s_kill formation
    γ_2::Real=0.68,                             # share of time investment in kid's s_kill formation
    τ_0::Real=0.10,                             # intercept of progressive income tax
    τ_1::Real=0.04,                             # slope of progressive income tax
    g::Real=0.03,                               # government welfare transfer
    q_A::Real=1.70,                             # adult consumption equivalance scale in data
    p_1::Real=0.33,                             # intercept of social security payment 
    p_0::Real=0.08,                             # slope of social security payment
    τ_s::Real=0.11,                             # payroll tax
    τ_k::Real=0.15,                             # interest income tax fpr retirees
    d_0::Real=0.02,                             # public subsidy in kid's education
    d_1::Real=0.09,                             # public subsidy in kid's education
    d_2::Real=0.10,                             # public subsidy in kid's education
    κ::Real=0.19,                               # college monetary cost
    κ_tilde::Real=2 / 3,                          # college time cost
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
    b::Real=0.81,                               # Ben-Porath human capital accumulation
    ζ::Real=3.11,                               # kid to adult human capital anchor
    ω_1::Real=0.56,                             # primary productivity of kid's investment
    ω_2::Real=0.30,                             # secondary productivity of kid's investment
    ψ_1::Real=0.23,                             # preference for kid to college with high school parent
    ψ_2::Real=0.24,                             # preference for kid to college with college parent
)
    """
    contruct an immutable object of all paramters
    """
    # adult consumption equivalance scale in model
    q = (1.0 + θ^(1.0 / χ))^χ

    # normal factor for utility
    normal_factor = (1.0 - β) * (1 - β^5 * θ) / (1.0 - β^9)

    # return values
    return (
        age_periods=age_periods,
        χ=χ,
        α=α,
        δ=δ,
        σ=σ,
        γ_0=γ_0,
        γ_1=γ_1,
        γ_2=γ_2,
        τ_0=τ_0,
        τ_1=τ_1,
        g=g,
        q=q,
        q_A=q_A,
        p_1=p_1,
        p_0=p_0,
        τ_s=τ_s,
        τ_k=τ_k,
        d_0=d_0,
        d_1=d_1,
        d_2=d_2,
        κ=κ,
        κ_tilde=κ_tilde,
        β=β,
        ν=ν,
        A=A,
        ϵ_μ=ϵ_μ,
        ϵ_σ=ϵ_σ,
        ϕ_0=ϕ_0,
        ϕ_1=ϕ_1,
        ϕ_2=ϕ_2,
        θ=θ,
        a_ρ=a_ρ,
        a_μ=a_μ,
        a_σ=a_σ,
        b=b,
        ζ=ζ,
        ω_1=ω_1,
        ω_2=ω_2,
        ψ_1=ψ_1,
        ψ_2=ψ_2,
        normal_factor=normal_factor,
    )
end