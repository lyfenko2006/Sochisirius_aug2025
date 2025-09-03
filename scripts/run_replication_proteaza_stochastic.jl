using Revise
using DrWatson
using DataFrames
using CSV
using DifferentialEquations
using JLD2

@quickactivate "Sochisirius_aug2025"
include(joinpath(@__DIR__, "../src/All_function.jl"))
using .All_function: jumps, vartojumps, jumptovars

V_mat_pi=[1.0, 0.516127832773321, 0.19430300396205602, 0.09185267311355345, 0.05179090287289477, 0.032763480953055577, 0.02241338560317059, 0.016216282455122053, 0.01223434324986414, 0.009534211199217702, 0.00762412134124905, 0.0062260564239387875, 0.005173677035519181, 0.0043627241508042295, 0.003725262985858844]

println("V_mat_pi=", V_mat_pi)

# parameters
pt = CSV.File(datadir("parameters.txt"), header=false, delim=' ', ignorerepeated=true)
pnames = String.(pt.Column1)
p = Float64.(pt.Column2)

# solve the model
t_final = 36.0
tspan = (0.0, t_final)

# initial condition
u0 = repeat([0.0], 25)
V0 = 10.0
u0[1] = V0

discrete_prob = DiscreteProblem(Int64.(u0), tspan, p)
jump_prob = JumpProblem(discrete_prob,
                            Direct(),
                            jumps...,
                            vartojumps_map=vartojumps,
                            jumptovars_map=jumptovars,
                            save_positions=(false,false))

solve(jump_prob, SSAStepper(), saveat=1.0)
n = 10000 + 1
vars_to_save = [6]
output_func(local_sol, i) = (local_sol[end][vars_to_save], false)
ens_prob = EnsembleProblem(jump_prob, output_func=output_func)

sim = solve(ens_prob, SSAStepper(),
                            EnsembleThreads(),
                            saveat=tspan[end],
                            trajectories=n)
 
p_infected = 1.0 - sum([u[1] for u in sim.u] .== 0) / n
coeff =  p_infected

println("Start calculate matrix")

t_final = 36
tspan = (0.0, t_final)
println(length(V_mat_pi))
matrix = zeros(length(V_mat_pi))

for (i, eps) in enumerate(V_mat_pi)
    local_p = copy(p)

    local_V0 = round(10*eps)

    local_u0 = repeat([0.0], 25)
    local_u0[1] = local_V0
    local_discrete_prob = DiscreteProblem(Int64.(local_u0), tspan, local_p)
    local_jump_prob = JumpProblem(local_discrete_prob,
                                        Direct(),
                                        jumps...,
                                        vartojumps_map=vartojumps,
                                        jumptovars_map=jumptovars,
                                        save_positions=(false,false))

    solve(local_jump_prob, SSAStepper(), saveat=1.0)

    local_n = 10000 + 1
    local_vars_to_save = [6]
    local_output_func(local_sol, i) = (local_sol[end][local_vars_to_save], false)
    local_ens_prob = EnsembleProblem(local_jump_prob, output_func=local_output_func)

    local_sim = solve(local_ens_prob, SSAStepper(),
                            EnsembleThreads(),
                            saveat=tspan[end],
                            trajectories=local_n)

    local_p_infected = 1.0 - sum([u[1] for u in local_sim.u] .== 0) / n
    matrix[i] = local_p_infected / coeff
    
end

println(size(matrix))
println(matrix)

data_dict = Dict("matrix" => matrix)
safesave(datadir("sims", "matrix_PI.h5"), data_dict)

println("Done!")