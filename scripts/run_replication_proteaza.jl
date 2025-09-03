using Revise
using DrWatson
using DataFrames
using CSV
using DifferentialEquations
using JLD2

@quickactivate "Sochisirius_aug2025"
include(joinpath(@__DIR__, "../src/All_function.jl"))
using .All_function: HIV_replication, inhibition, apply_inhibitors_1!
# Данные для DRV (PI)
IC50_DRV = 0.0264 * 547.66 / 1000

# inhibitors
inhibitors = DataFrame(
    tag = ["DRV"],
    IC50 = [IC50_DRV],
    m = [2.14],
    target = ["k_mat"] 
)

pt = CSV.File(datadir("parameters.txt"), header=false, delim=' ', ignorerepeated=true)
pnames = String.(pt.Column1)
p = Float64.(pt.Column2)
lb = Float64.(pt.Column3)
ub = Float64.(pt.Column4)
@assert all(lb .<= p .<= ub)

# initial condition
u0 = repeat([0.0],26)
V0 = 10.0
u0[1] = V0
t_final = 36
tspan = (0.0, t_final)

println("Start calculate coeff")
prob = ODEProblem(HIV_replication,u0,tspan,p)
sol = solve(prob,RadauIIA5(),reltol=1e-14,abstol=1e-14)
var12 = [u[12] for u in sol.u]
coeff = (var12[length(var12)] + var12[1])/2 * t_final

u24 = [u[24] for u in sol.u]
V_mat = u24[length(u24)]
println("coeff=", coeff, " V_mat=", V_mat)

println("Start calculate matrix")

V_mat_pi = []
for (i, conc1) in enumerate(0:0.005:0.07)
    local_p = copy(p)
    apply_inhibitors_1!(local_p, inhibitors, conc1)
            
    local_prob = ODEProblem(HIV_replication,u0,tspan,local_p)
    local_sol = solve(local_prob,RadauIIA5(),reltol=1e-14,abstol=1e-14)
    local_u24 = [u[24] for u in local_sol.u]
    temp = local_u24[length(local_u24)]
    push!(V_mat_pi, temp/V_mat)
end
 
println("V_mat_pi=", V_mat_pi)

println("Done!")