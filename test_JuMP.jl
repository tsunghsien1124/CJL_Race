using JuMP
import Ipopt

#==============#
# Syntax notes #
#==============#
model = Model();
@variable(model, x[1:2]);
@variable(model, y);
c = [1, 2];
# incorrect way
@NLobjective(model, Min, c' * x + 3y) # Unexpected array [1 2] in nonlinear expression. Nonlinear expressions may contain only scalar expressions.
# using sum function
@NLobjective(model, Min, sum(c[i] * x[i] for i = 1:2) + 3y)
# using @expression
@expression(model, expr, c' * x)
@NLobjective(model, Min, expr + 3y)

#===========#
# Splatting #
#===========#
model = Model();
@variable(model, x[1:3]);
@NLconstraint(model, *(x...) <= 1.0)
# incorrect way
@NLconstraint(model, *((x / 2)...) <= 0.0) # ERROR: Unsupported use of the splatting operator. JuMP supports splatting only symbols. For example, `x...` is ok, but `(x + 1)...`, `[x; y]...` and `g(f(y)...)` are not.

#========================#
# User-defined Functions #
#========================#
f(x::Float64) = 2 * x  # This will not work.
f(x::Real)    = 2 * x  # This is good.
f(x)          = 2 * x  # This is also good.

function bad_f(x...)
    y = zeros(length(x))  # This constructs an array of `Float64`!
    for i = 1:length(x)
        y[i] = x[i]^i
    end
    return sum(y)
end

function good_f(x::T...) where {T<:Real}
    y = zeros(T, length(x))  # Construct an array of type `T` instead!
    for i = 1:length(x)
        y[i] = x[i]^i
    end
    return sum(y)
end

using Optim
f(x) = (x[1] - 2.0)^2 + (x[2] - 1.0)^2 + 1.0
res = optimize(f, [0.0, 0.0], [3.0, 3.0], [2.9, 2.9], Fminbox(BFGS()))
println(res.minimum)  # Prints 10.250625
println(res.minimizer) # Prints [4.375, 2.9]