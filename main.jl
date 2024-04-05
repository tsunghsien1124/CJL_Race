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

#======#
# Test #
#======#
parameters = parameters_function();
prices = prices_function(parameters);
variables = variables_function(prices, parameters);

age_10_function!(variables, prices, parameters);
age_9_function!(variables, prices, parameters);
age_8_function!(variables, prices, parameters);
age_7_function!(variables, prices, parameters);
age_6_function!(variables, prices, parameters);
