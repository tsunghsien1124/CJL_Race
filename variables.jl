mutable struct Mutable_Variables
    """
    construct a type for mutable variables
    """
    V_4::Array{Float64,4}
    policy_s_5::Array{Float64,4}
    policy_n_4::Array{Float64,4}
    V_5::Array{Float64,5}
    policy_s_6::Array{Float64,5}
    policy_n_5::Array{Float64,5}
    policy_l_0::Array{Float64,5}
    V_6::Array{Float64,6}
    policy_s_7::Array{Float64,6}
    policy_n_6::Array{Float64,6}
    policy_l_1::Array{Float64,6}
    V_7::Array{Float64,6}
    policy_s_8::Array{Float64,6}
    policy_n_7::Array{Float64,6}
    policy_l_2::Array{Float64,6}
    V_8::Array{Float64,6}
    W_8::Array{Float64,7}
    policy_s_9::Array{Float64,7}
    policy_n_8::Array{Float64,7}
    policy_n_3::Array{Float64,7}
    policy_c_3::Array{Float64,6}
    V_9::Array{Float64,7}
    policy_s_10::Array{Float64,7}
    policy_n_9::Array{Float64,7}
    policy_s_4::Array{Float64,7}
    V_10::Array{Float64,3}
end

function variables_function(parameters::NamedTuple)
    """
    construct a mutable object containing endogenous variables
    """
    @unpack c_size, a_size, h_size, s_size = parameters

    # j = 4 (independece)
    V_4 = zeros(c_size, a_size, h_size, s_size)
    policy_s_5 = zeros(c_size, a_size, h_size, s_size)
    policy_n_4 = zeros(c_size, a_size, h_size, s_size)

    # j = 5 (child investment)
    V_5 = zeros(a_size, c_size, a_size, h_size, s_size)
    policy_s_6 = zeros(a_size, c_size, a_size, h_size, s_size)
    policy_n_5 = zeros(a_size, c_size, a_size, h_size, s_size)
    policy_l_0 = zeros(a_size, c_size, a_size, h_size, s_size)

    # j = 6
    V_6 = zeros(a_size, h_size, c_size, a_size, h_size, s_size)
    policy_s_7 = zeros(a_size, h_size, c_size, a_size, h_size, s_size)
    policy_n_6 = zeros(a_size, h_size, c_size, a_size, h_size, s_size)
    policy_l_1 = zeros(a_size, h_size, c_size, a_size, h_size, s_size)

    # j = 7
    V_7 = zeros(a_size, h_size, c_size, a_size, h_size, s_size)
    policy_s_8 = zeros(a_size, h_size, c_size, a_size, h_size, s_size)
    policy_n_7 = zeros(a_size, h_size, c_size, a_size, h_size, s_size)
    policy_l_2 = zeros(a_size, h_size, c_size, a_size, h_size, s_size)

    # j = 8 (college decision)
    V_8 = zeros(a_size, h_size, c_size, a_size, h_size, s_size)
    W_8 = zeros(c_size, a_size, h_size, c_size, a_size, h_size, s_size)
    policy_s_9 = zeros(c_size, a_size, h_size, c_size, a_size, h_size, s_size)
    policy_n_8 = zeros(c_size, a_size, h_size, c_size, a_size, h_size, s_size)
    policy_n_3 = zeros(c_size, a_size, h_size, c_size, a_size, h_size, s_size)
    policy_c_3 = zeros(a_size, h_size, c_size, a_size, h_size, s_size)

    # j = 9 (inter vivos)
    V_9 = zeros(c_size, a_size, h_size, c_size, a_size, h_size, s_size)
    policy_s_10 = zeros(c_size, a_size, h_size, c_size, a_size, h_size, s_size)
    policy_n_9 = zeros(c_size, a_size, h_size, c_size, a_size, h_size, s_size)
    policy_s_4 = zeros(c_size, a_size, h_size, c_size, a_size, h_size, s_size)

    # j = 10 (deterministic)
    V_10 = zeros(c_size, h_size, s_size)

    # return all variable placeholders
    variables = Mutable_Variables(
        V_4,
        policy_s_5,
        policy_n_4,
        V_5,
        policy_s_6,
        policy_n_5,
        policy_l_0,
        V_6,
        policy_s_7,
        policy_n_6,
        policy_l_1,
        V_7,
        policy_s_8,
        policy_n_7,
        policy_l_2,
        V_8,
        W_8,
        policy_s_9,
        policy_n_8,
        policy_n_3,
        policy_c_3,
        V_9,
        policy_s_10,
        policy_n_9,
        policy_s_4,
        V_10,
    )
    return variables
end

mutable struct Mutable_Prices
    """
    construct a type for mutable prices (to be solved in general equlibrium)
    """
    W_S::Vector{Float64}
end