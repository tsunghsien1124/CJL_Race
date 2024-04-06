#=================#
# Import Packages #
#=================#
using Distributions: Normal, quantile
using QuantEcon: rouwenhorst, stationary_distributions, MarkovChain
using Parameters: @unpack
using DelimitedFiles
# using JuMP
# import Ipopt
using Optim
using ProgressMeter
using BenchmarkTools
using Plots
using JLD2: @save, @load

#==================#
# Import Functions #
#==================#
include("functions.jl")
include("parameters.jl")
include("variables.jl")
include("age_10.jl")
include("age_9.jl")
include("age_8.jl")
include("age_7.jl")
include("age_6.jl")
include("age_5.jl")
include("age_4.jl")

#================#
# Initialization #
#================#
parameters = parameters_function();
prices = prices_function(parameters);
variables = variables_function(prices, parameters);
load_initial_V_4 = true

#==========================#
# Value Function Iteration #
#==========================#
age_10_function!(variables, prices, parameters); # age 10 is deterministic
if load_initial_V_4 == true
    @load "V_4_temp.jld2" V_4_temp
    copyto!(variables.V_4, V_4_temp)
end
V_4_temp = similar(variables.V_4);
crit = 1E-6;
diff = Inf;
while diff > crit
    copyto!(V_4_temp, variables.V_4);
    age_9_function!(variables, prices, parameters);
    age_8_function!(variables, prices, parameters);
    age_7_function!(variables, prices, parameters);
    age_6_function!(variables, prices, parameters);
    age_5_function!(variables, prices, parameters);
    age_4_function!(variables, prices, parameters);
    diff = maximum(abs.(V_4_temp .- variables.V_4));
    println(diff)
end
@save "V_4_temp.jld2" V_4_temp

#=============#
# Check Plots #
#=============#
# age 10
h_ind_age_10 = parameters.h_size
plot(parameters.s_grid, variables.V_10[:, h_ind_age_10, 1], label="high school")
plot!(parameters.s_grid, variables.V_10[:, h_ind_age_10, 2], label="college")