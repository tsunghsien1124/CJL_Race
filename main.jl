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
variables = variables_function(parameters);
prices = prices_function(parameters);

age_10_function!(variables, prices, parameters)