using Revise
using DrWatson
using DataFrames
using CSV
using DifferentialEquations
using JLD2
using Statistics
using Base.Threads
using SparseArrays

@quickactivate "Sochisirius_aug2025"
include(joinpath(@__DIR__, "../src/All_function.jl"))
using .All_function: apply_inhibitors_3!, jumps, vartojumps, jumptovars

Threads.nthreads() = max(4, Sys.CPU_THREADS) 
println("Доступно потоков: ", nthreads())

# V_mat_pi=[1.0, 0.516127832773321, 0.19430300396205602, 0.09185267311355345, 0.05179090287289477, 0.032763480953055577, 0.02241338560317059, 0.016216282455122053, 0.01223434324986414, 0.009534211199217702, 0.00762412134124905, 0.0062260564239387875, 0.005173677035519181, 0.0043627241508042295, 0.003725262985858844]

# println("V_mat_pi=", V_mat_pi)

# inhibitors
inhibitors = DataFrame(
    tag = ["DRV"],
    IC50 = [0.0146],
    m = [2.14],
    target = ["k_mat"] 
)

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

# parameters
pt = CSV.File(datadir("parameters.txt"), header=false, delim=' ', ignorerepeated=true)
pnames = String.(pt.Column1)
p = Float64.(pt.Column2)

# solve the model
t_final = 36.0
tspan = (0.0, t_final)

# initial condition
# u0 = repeat([0.0], 25)
# V0 = 4.0
# u0[1] = V0

println("Start calculate matrix")
V_mat_pi = [1.0, 0.05179090287289477, 0.01223434324986414, 0.005173677035519181, 
            0.002801891631620362, 0.00173988377236717, 0.0011784561020117895, 
            0.0008475967122302999, 0.0006370548450153384, 0.0004951893284051425, 0.0003952693153442152]

println("V_mat_pi=", V_mat_pi)
t_final = 36
tspan = (0.0, t_final)
count_calc = 100

NRTI_conc_range = 0.0:0.1:1  # значения для 3TC и AZT
n_NRTI = length(NRTI_conc_range)

matrix = zeros(count_calc, 11, n_NRTI, n_NRTI)

Threads.@threads for j in 1:count_calc
    println("Thread $(Threads.threadid()) is processing column j=$j")
    local_matrix = zeros(11, n_NRTI, n_NRTI)  # Локальный массив для каждого потока

    for (i, eps) in enumerate(V_mat_pi)
        # Начальные условия считаем один раз:
        local_u0 = repeat([0.0], 26)
        v = rand(10)
        local_u0[1] = sum(x -> x < eps, v)

        for (idx2, conc2) in enumerate(NRTI_conc_range)
            for (idx3, conc3) in enumerate(NRTI_conc_range)
                
                local_p = copy(p)
                apply_inhibitors_3!(local_p, inhibitors, 1.0, conc2, conc3)  # DRV=0.0, NRTIs varying

                discrete_prob = DiscreteProblem(Int64.(local_u0), tspan, local_p)
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
                            EnsembleSerial(),
                            saveat=tspan[end],
                            trajectories=n)

                p_infected = 1.0 - sum([u[1] for u in sim.u] .== 0) / n
                local_matrix[i, idx2, idx3] = p_infected
            end
        end
    end
   
    # Синхронизированная запись в общий массив
    @inbounds for i in 1:11, idx2 in 1:n_NRTI, idx3 in 1:n_NRTI
        matrix[j, i, idx2, idx3] = local_matrix[i, idx2, idx3]
    end
end

println(size(matrix))

data_dict = Dict("matrix" => matrix)
safesave(datadir("sims", "matrix3d_$count_calc.h5"), data_dict)

println("Done!")