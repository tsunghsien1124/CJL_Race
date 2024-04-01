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

#==================#
# Import Functions #
#==================#
include("functions.jl")
include("parameters.jl")
include("variables.jl")
include("age_10.jl")
include("age_9.jl")

#======#
# Test #
#======#
parameters = parameters_function();
prices = prices_function(parameters);
variables = variables_function(prices, parameters);

age_10_function!(variables, prices, parameters);
age_9_function!(variables, prices, parameters);