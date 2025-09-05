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
using .All_function: jumps, vartojumps, jumptovars

Threads.nthreads() = max(4, Sys.CPU_THREADS) 
println("Доступно потоков: ", nthreads())

# V_mat_pi=[1.0, 0.516127832773321, 0.19430300396205602, 0.09185267311355345, 0.05179090287289477, 0.032763480953055577, 0.02241338560317059, 0.016216282455122053, 0.01223434324986414, 0.009534211199217702, 0.00762412134124905, 0.0062260564239387875, 0.005173677035519181, 0.0043627241508042295, 0.003725262985858844]

# println("V_mat_pi=", V_mat_pi)

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

println("Start calculate matrix")

V_mat_pi = [1.0, 0.10272495445905991, 0.025290524097613155, 
        0.010775890824127686, 0.005850741311279047, 0.0036372714423588883, 
        0.0024650813357549586, 0.001773624961794668, 0.001333361192681484, 
        0.0010365938286168179, 0.0008275174942637705]

println("V_mat_pi=", V_mat_pi)
t_final = 36
tspan = (0.0, t_final)
count_calc = 1000
matrix = zeros(11, count_calc)

Threads.@threads for j in 1:count_calc
    # println("Thread $(Threads.threadid()) is processing column j=$j")
    local_matrix = zeros(11)  # Локальный массив для каждого потока
    for (i, eps) in enumerate(V_mat_pi)
        local_p = copy(p)
        local_u0 = repeat([0.0], 26)
        v = rand(400)
        local_u0[1] = sum(x -> x < eps, v)
        # println("j=$j, i=$i, V0=", local_u0[1])

        discrete_prob = DiscreteProblem(Int64.(local_u0), tspan, local_p)
        jump_prob = JumpProblem(discrete_prob,
                                Direct(),
                                jumps...,
                                vartojumps_map=vartojumps,
                                jumptovars_map=jumptovars,
                                save_positions=(false,false))
        @timed local_sol = solve(jump_prob, SSAStepper(), saveat=1.0)

        n = 10000 + 1
        vars_to_save = [6]
        output_func(local_sol, i) = (local_sol[end][vars_to_save], false)
        ens_prob = EnsembleProblem(jump_prob, output_func=output_func)

        @timed sim = solve(ens_prob, SSAStepper(),
                          EnsembleSerial(),
                          saveat=tspan[end],
                          trajectories=n)

        p_infected = 1.0 - sum([u[1] for u in sim.u] .== 0) / n
        # display(p_infected)
        local_matrix[i] = p_infected
    end

    # Синхронизированная запись в общий массив
    @inbounds for i in 1:11
        matrix[i, j] = local_matrix[i]
    end
end

using SparseArrays, Statistics

row_means = [mean(row) for row in eachrow(matrix)]


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

print_matrix(matrix)
println(row_means)

# h5open("data/matrix_$count_calc.h5", "w") do file
#     write(file, "matrix", matrix)
# end
data_dict = Dict("matrix" => matrix)
safesave(datadir("sims", "matrix_$count_calc.h5"), data_dict)
# h5open("data/mean_values_$count_calc.h5", "w") do file
#     write(file, "mean_values", row_means)
# end
data_dict = Dict("mean_values" => row_means)
safesave(datadir("sims", "mean_values_$count_calc.h5"), data_dict)
println("Done!")