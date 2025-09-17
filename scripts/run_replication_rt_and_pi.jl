using Revise
using DrWatson
using DataFrames
using CSV
using DifferentialEquations
using JLD2

@quickactivate "Sochisirius_aug2025"
include(joinpath(@__DIR__, "../src/All_function.jl"))
using .All_function: HIV_replication, inhibition, apply_inhibitors!

# Данные для 3TC и AZT (NRTI)
IC50_3TC = 0.00363258 # µM
IC50_AZT = 0.023487990000000004  # µM
# Данные для DRV (PI)
IC50_DRV = 0.0264 * 547.66 / 1000 / 7

# inhibitors
inhibitors = DataFrame(
    tag = ["DRV", "3TC", "AZT"],
    IC50 = [IC50_DRV, IC50_3TC, IC50_AZT],
    m = [2.14, 0.97, 0.62],
    target = ["k_mat", "k_RT", "k_RT"] 
)

pt = CSV.File(datadir("parameters.txt"), header=false, delim=' ', ignorerepeated=true)
pnames = String.(pt.Column1)
p = Float64.(pt.Column2)
lb = Float64.(pt.Column3)
ub = Float64.(pt.Column4)
@assert all(lb .<= p .<= ub)

# initial condition
u0 = repeat([0.0],27)
V0 = 4.0
u0[1] = V0
t_final = 36.0
tspan = (0.0, t_final)

println("Start calculate coeff")
prob = ODEProblem(HIV_replication,u0,tspan,p)
sol = solve(prob,RadauIIA5(),reltol=1e-14,abstol=1e-14)
var27 = [u[27] for u in sol.u]
coeff = (var27[length(var27)] + var27[1])/2 * t_final 
println("coeff=", coeff)

println("Start calculate matrix")
matrix = zeros(15, 11, 11)

for (i, conc1) in enumerate(0:0.005:0.07)
    for (j, conc2) in enumerate(0:0.1:1)
        for (k, conc3) in enumerate(0:0.1:1)
            local_u0 = repeat([0.0], 27)
            local_u0[1] = V0
            local_p = copy(p)
            concentrations = [conc1, conc2, conc3]
            apply_inhibitors!(local_p, inhibitors, concentrations)
            
            local_prob = ODEProblem(HIV_replication,local_u0,tspan,local_p)
            local_sol = solve(local_prob,RadauIIA5(),reltol=1e-14,abstol=1e-14)

            local_u27 = [u[27] for u in local_sol.u]
            matrix[i, j, k] = local_u27[length(local_u27)] / coeff
        end
    end
end

data_dict = Dict("matrix" => matrix)
safesave(datadir("sims", "matrix_u2.h5"), data_dict)
println("Done!")