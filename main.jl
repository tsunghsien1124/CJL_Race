#=================#
# Import Packages #
#=================#
using Distributions: Normal, quantile
using QuantEcon: rouwenhorst, stationary_distributions, MarkovChain, Categorical
using Parameters: @unpack
# using DelimitedFiles
# using JuMP
# import Ipopt
using Optim
using ProgressMeter
using BenchmarkTools
using Plots
using JLD2: @save, @load
using Random

#==================#
# Import Functions #
#==================#
include("parameters.jl")
include("variables.jl")
include("grids.jl")
include("functions.jl")
include("age_10.jl")
include("age_9.jl")
include("age_8.jl")
include("age_7.jl")
include("age_6.jl")
include("age_5.jl")
include("age_4.jl")
include("Simulation.jl")

#================#
# Initialization #
#================#
parameters = parameters_function();
prices = prices_function(parameters);
grids = grids_function(parameters, prices);
variables = variables_function(prices, parameters);

#==========================#
# Value Function Iteration #
#==========================#
function VFI!(variables::Mutable_Variables, prices::Mutable_Prices, parameters::NamedTuple, grids::NamedTuple; crit = 1E-6, diff = Inf, load_initial_V_4 = true)
    """
    conduct backward value function iterations
    """
    age_10_function!(variables, prices, parameters, grids); # age 10 is deterministic
    if load_initial_V_4 == true
        @load "V_4_temp.jld2" V_4_temp
        copyto!(variables.V_4, V_4_temp)
    end
    V_4_temp = similar(variables.V_4);
    while diff > crit
        copyto!(V_4_temp, variables.V_4);
        age_9_function!(variables, prices, parameters, grids);
        age_8_function!(variables, prices, parameters, grids);
        age_7_function!(variables, prices, parameters, grids);
        age_6_function!(variables, prices, parameters, grids);
        age_5_function!(variables, prices, parameters, grids);
        age_4_function!(variables, prices, parameters, grids);
        diff = maximum(abs.(V_4_temp .- variables.V_4));
        println(diff)
    end
    @save "V_4_temp.jld2" V_4_temp
end
VFI!(variables, prices, parameters, grids)

#============#
# Simulation #
#============#


#=============#
# Check Plots #
#=============#
# age 10
h_ind_age_10 = grids.h_size
plot(grids.s_grid, variables.V_10[:, h_ind_age_10, 1], label="high school")
plot!(grids.s_grid, variables.V_10[:, h_ind_age_10, 2], label="college")