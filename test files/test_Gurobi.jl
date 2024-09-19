using JuMP
import Gurobi
import Interpolations
import PiecewiseLinearOpt

I = [1.0, 1.2, 1.3]
J = [1.0, 1.2, 1.3]
V = rand(3, 3)
interp1 = Interpolations.interpolate(
    (I, J),
    V,
    Interpolations.Gridded(Interpolations.Linear()),
)
W = rand(3, 3)
interp2 = Interpolations.interpolate(
    (I, J),
    W,
    Interpolations.Gridded(Interpolations.Linear()),
)
model = Model(Gurobi.Optimizer)
@variable(model, 0.0 <= c1 <= 1.0)
@variable(model, 0.0 <= c2 <= 1.0)
@variable(model, 0.0 <= c3 <= 1.0)
z1 = PiecewiseLinearOpt.piecewiselinear(model, c1, c2, I, J, (u,v) -> interp1[u, v])
z2 = PiecewiseLinearOpt.piecewiselinear(model, c1, c3, I, J, (u,v) -> interp2[u, v])
@objective(model, Min,  1.5c1 + 1.3c2 + z1 + z2)
optimize!(model)