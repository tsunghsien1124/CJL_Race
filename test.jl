using Optimization
using OptimizationNLopt
# rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
# x0 = zeros(2)
# p = [1.0, 100.0]
# f = OptimizationFunction(rosenbrock)
# prob = Optimization.OptimizationProblem(f, x0, p, lb = [-1.0, -1.0], ub = [1.0, 1.0])
# sol = solve(prob, NLopt.LN_NELDERMEAD())

using Interpolations
import Random
Random.seed!(1234)

function test()
    # Sample data points for the interpolation
    α = 2.0
    data_1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    data_2 = [1.0, 1.5, 2.0, 2.1, 1.4, 0.2]
    itp = linear_interpolation(data_1, data_2, extrapolation_bc=Line()) 
    f_itp(x) = -itp(x[1])
    function f_obj(x, p)
        if p[1] >= 0.8
            return (x[1]-5.2)^2 + log(max(p[1], p[2]))
        else
            return (x[1]-2.8)^2
        end
    end
    f_OBJ(x, p) = -itp(x[1]) + f_obj(x, p) + α
    p_vec = rand(2,4)
    x0 = zeros(1)
    results_x = zeros(size(p_vec, 2))
    results_obj = zeros(size(p_vec, 2))
    f = OptimizationFunction(f_OBJ)
    for p_i in 1:size(p_vec, 2)
        itp.itp.coefs .= data_2 * (1.0 + p_i/100)
        p = p_vec[:, p_i]
        prob = Optimization.OptimizationProblem(f, x0, p, lb = [-200], ub = [14])
        sol = solve(prob, NLopt.LN_COBYLA())
        results_x[p_i] = sol.u[1]
        results_obj[p_i] = sol.objective
    end
    return results_x, results_obj
end