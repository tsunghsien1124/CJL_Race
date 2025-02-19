function grids_function(parameters::NamedTuple, prices::Mutable_Prices;
    #-------------------------#
    # numerical specification #
    #-------------------------#
    age_min::Int64=4,                         # min model age
    age_max::Int64=12,                        # max model age
    c_size::Int64=2,                          # number of academic degree 
    ϵ_size::Int64=2,                          # number of market luck shock
    a_size::Int64=3,                          # number of learning ability shock
    h_size::Int64=6,                          # number of human capital
    h_power::Int64=2,                         # grid power of human capital
    h_frac::Float64=20.0,                     # upper bound scale of human capital
    h_k_frac::Float64=0.25,                   # upper bound scale of kid's human capital
    s_size::Int64=12,                         # number of asset choice
    s_power::Int64=2,                         # grid power of asset choice
    s_frac::Float64=6.0,                      # upper bound scale of asset choices
    s_k_frac::Float64=1.25,                   # upper bound scale of kid's asset choices
)
    """
    contruct an immutable object of all grids
    """
    @unpack α, g, ϵ_μ, ϵ_σ, a_ρ, a_σ, a_μ = parameters

    # age
    age_size = length(age_min:age_max)

    # market luck shock
    ϵ_grid = ϵ_grid_function(ϵ_size, ϵ_μ, ϵ_σ, c_size)

    # learning ability shock
    a_Γ, a_grid, a_birth = a_grid_function(a_size, a_ρ, a_σ, a_μ)

    # human capital (parent and kid)
    h_min = 0.0
    h_max = prices.e_bar / prices.W * h_frac
    h_step, h_grid = h_grid_function(h_size, h_min, h_max, h_power)

    h_k_max = h_max * h_k_frac
    h_k_step, h_k_grid = h_grid_function(h_size, h_min, h_k_max, h_power)

    # saving or borrowing (parent and kid)
    s_min = -g / (1.0 + prices.r)
    s_max = prices.e_bar * s_frac
    s_step, s_grid = h_grid_function(s_size, s_min, s_max, s_power)

    s_k_min = 0.0
    s_k_max = s_max * s_k_frac
    s_k_step, s_k_grid = h_grid_function(s_size, s_k_min, s_k_max, s_power)

    # labor choice
    # n_step = 0.1
    # n_grid = collect(0.0:n_step:1.0)
    # n_size = length(n_grid)

    # loop index
    loop_age_10 = collect(Iterators.product(1:c_size, 1:a_size, 1:h_size, 1:s_size, 1:c_size, 1:a_size, 1:h_size))

    # return values
    return (
        age_min=age_min,
        age_max=age_max,
        age_size=age_size,
        c_size=c_size,
        ϵ_size=ϵ_size,
        ϵ_grid=ϵ_grid,
        a_size=a_size,
        a_Γ=a_Γ,
        a_grid=a_grid,
        a_birth=a_birth,
        h_size=h_size,
        h_power=h_power,
        h_frac=h_frac,
        h_k_frac=h_k_frac,
        h_min=h_min,
        h_max=h_max,
        h_step=h_step,
        h_grid=h_grid,
        h_k_max=h_k_max,
        h_k_step=h_k_step,
        h_k_grid=h_k_grid,
        s_size=s_size,
        s_power=s_power,
        s_frac=s_frac,
        s_k_frac=s_k_frac,
        s_min=s_min,
        s_max=s_max,
        s_step=s_step,
        s_grid=s_grid,
        s_k_min=s_k_min,
        s_k_max=s_k_max,
        s_k_step=s_k_step,
        s_k_grid=s_k_grid,
        # n_step = n_step,
        # n_grid = n_grid,
        # n_size = n_size,
        loop_age_10 = loop_age_10,
    )
end