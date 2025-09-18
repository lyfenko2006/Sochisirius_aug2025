using Revise
using DrWatson
using DifferentialEquations
using CSV, DataFrames
using HDF5 

@quickactivate "Sochisirius_aug2025"
include(joinpath(@__DIR__, "../src/All_function.jl"))
using .All_function: inhibition, apply_inhibitors!, jumps, vartojumps, jumptovars

# Данные для 3TC и AZT (NRTI)
IC50_3TC = 0.003633 # µM
IC50_AZT = 0.0235 # µM

# inhibitors
inhibitors = DataFrame(
    tag = ["3TC", "AZT"],
    IC50 = [IC50_3TC, IC50_AZT],
    m = [0.97, 0.62],
    target = ["k_RT", "k_RT"] 
)

# parameters
pt = CSV.File(datadir("parameters.txt"), header=false, delim=' ', ignorerepeated=true)
pnames = String.(pt.Column1)
p = Float64.(pt.Column2)

# solve the model
t_final = 36.0
tspan = (0.0, t_final)

# initial condition
u0 = repeat([0.0], 25)
V0 = 4.0
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

NRTI_conc_range = 0.0:0.1:1 
n_NRTI = length(NRTI_conc_range)
vec_3TC = zeros(1000, n_NRTI)
vec_AZT = zeros(1000, n_NRTI)
for i in 1:1000
    for (idx2, conc2) in enumerate(NRTI_conc_range)
    conc3 = 0
    local_p = copy(p)
    concentrations = [conc2, conc3]
    apply_inhibitors!(local_p, inhibitors, concentrations) 
    local_discrete_prob = DiscreteProblem(Int64.(u0), tspan, local_p)
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
    vec_3TC[i, idx2] = local_p_infected / coeff
    end

    for (idx3, conc3) in enumerate(NRTI_conc_range)
        conc2 = 0
        local_p = copy(p)
        concentrations = [conc2, conc3]
        apply_inhibitors!(local_p, inhibitors, concentrations) 
        local_discrete_prob = DiscreteProblem(Int64.(u0), tspan, local_p)
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
        vec_AZT[i, idx3] = local_p_infected / coeff
    end
end

println(size(vec_AZT))

using Printf

function print_matrix(matrix, conc_range=0:0.1:1)
    # Шапка с концентрациями
    @printf("%6s", "conc1\\conc2")
    for conc in conc_range
        @printf("%8.2f", conc)
    end
    println()

    # Тело матрицы
    for (i, conc1) in enumerate(conc_range)
        @printf("%6.2f", conc1) 
        for j in 1:size(matrix, 2)
            @printf("%8.4f", matrix[i, j])
        end
        println()
    end
end

function print_vec(vec, conc_range=0:0.1:1)
    # Шапка с концентрациями
    # @printf("%6s", "conc1\\conc2")
    for conc in conc_range
        @printf("%8.2f", conc)
    end
    println()

    # Тело матрицы
    for i in 1:size(vec, 1)
        # @printf("%6.2f", conc1) 
        for j in 1:size(vec, 2)
            @printf("%8.4f", vec[i, j])
        end
        println()
    end
end
vec_3TC = sort(vec_3TC, dims=1)
vec_AZT = sort(vec_AZT, dims=1)
print_vec(vec_3TC)
print_vec(vec_AZT)

data_dict = Dict("matrix" => vec_3TC)
safesave(datadir("sims", "vec_3TC_u1.h5"), data_dict)
data_dict = Dict("matrix" => vec_AZT)
safesave(datadir("sims", "vec_AZT_u1.h5"), data_dict)

println("Done!")