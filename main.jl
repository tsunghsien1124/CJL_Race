#=================#
# Import Packages #
#=================#
using Distributions: Normal, quantile
using QuantEcon: rouwenhorst, stationary_distributions, MarkovChain
using Parameters: @unpack
using DelimitedFiles

#==================#
# Import Functions #
#==================#
include("functions.jl")
include("parameters.jl")
include("variables.jl")
include("old_parent.jl")

#======#
# Test #
#======#
parameters = parameters_function();
variabels = variables_function(parameters);