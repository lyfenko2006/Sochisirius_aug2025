using DifferentialEquations
using CSV, DataFrames
using HDF5 

jumps = []

# mapping the variables to the rates that depend on them
vartojumps = Vector{Vector{Int64}}()

# mapping the jump events to the variables which are changed in them
jumptovars = Vector{Vector{Int64}}()
 
push!(vartojumps, [1,2]) # u[1] (rate1 and rate2 depend on u[1])
push!(vartojumps, [3,4]) # u[2] (rates 3 and 4 depend on u[2])
# TODO: ... push values to vartojumps for the rates 5-11 ...
#
# ...
push!(vartojumps, [3,5,6]) # u[3] (rates 3, 5, 6 depend on u[3])
push!(vartojumps, [5,7,8]) # u[4] (rates 5, 7 and 8 depend on u[4])
push!(vartojumps, [7,9,10]) # u[5] (rates 7, 9 and 10 depend on u[5])
push!(vartojumps, [9,11]) # u[6] (rates 9 and 11 depend on u[6])

#    #dV_free
#    du[1] = -k_bound*V_free - d*V_free
function rate1(u,p,t) 
    #k_bound*V_free
    p[1]*u[1]
end
function affect1!(integrator)
    integrator.u[1] -= 1
    integrator.u[2] += 1
    nothing
end
push!(jumps, ConstantRateJump(rate1, affect1!))
push!(jumptovars, [1,2]) # (affect1 changes u[1] and u[2])

function rate2(u,p,t)
    #d*V_free
    p[50]*u[1]
end
function affect2!(integrator)
    integrator.u[1] -= 1
end
push!(jumps, ConstantRateJump(rate2,affect2!))
push!(jumptovars, [1]) # (affect2 changes u[1])

#    #dV_bound
#    du[2] = k_bound*V_free - (k_fuse+d_bound)*V_bound
function rate3(u,p,t)
    #k_fuse*V_bound
    p[3]*u[2]
end
function affect3!(integrator)
    integrator.u[2] -= 1
    integrator.u[3] += 1
end
push!(jumps, ConstantRateJump(rate3,affect3!))
push!(jumptovars, [2,3]) # (affect3 changes u[2] and u[3])

function rate4(u,p,t)
    #d_bound*V_bound
    p[2]*u[2]
end
function affect4!(integrator)
    integrator.u[2] -= 1
end
push!(jumps, ConstantRateJump(rate4,affect4!))
push!(jumptovars, [2]) # (affect4 changes u[2])

# TODO: add rates 5 to 11
#
# ...
function rate5(u,p,t)
    #k_RT*RNA_cor
    p[4]*u[3]
end
function affect5!(integrator)
    integrator.u[3] -= 1
    integrator.u[4] += 1
end
push!(jumps, ConstantRateJump(rate5,affect5!))
push!(jumptovars, [3,4]) # (affect5 changes u[3] and u[4])

function rate6(u,p,t)
    #d_RNA_cor*RNA_cor
    p[5]*u[3]
end
function affect6!(integrator)
    integrator.u[3] -= 1
end
push!(jumps, ConstantRateJump(rate6,affect6!))
push!(jumptovars, [3]) # (affect6 changes u[3])

function rate7(u,p,t)
    #k_DNA_t*DNA_cor
    p[7]*u[4]
end
function affect7!(integrator)
    integrator.u[4] -= 1
    integrator.u[5] += 1
end
push!(jumps, ConstantRateJump(rate7,affect7!))
push!(jumptovars, [4,5]) # (affect7 changes u[4] and u[5])

function rate8(u,p,t)
    #d_DNA_cor * DNA_cor
    p[6]*u[4]
end
function affect8!(integrator)
    integrator.u[4] -= 1
end
push!(jumps, ConstantRateJump(rate8,affect8!))
push!(jumptovars, [4]) # (affect8 changes u[4])

function rate9(u,p,t)
    #k_int * DNA_nuc
    p[9]*u[5]
end
function affect9!(integrator)
    integrator.u[5] -= 1
    integrator.u[6] += 1
end
push!(jumps, ConstantRateJump(rate9,affect9!))
push!(jumptovars, [5,6]) # (affect9 changes u[5] and u[6])

function rate10(u,p,t)
    #dDNA_nuc * DNA_nuc
    p[8]*u[5]
end
function affect10!(integrator)
    integrator.u[5] -= 1
end
push!(jumps, ConstantRateJump(rate10,affect10!))
push!(jumptovars, [5]) # (affect10 changes u[5])

function rate11(u,p,t)
    #dDNA_int * DNA_int
    p[10]*u[6]
end
function affect11!(integrator)
    integrator.u[6] -= 1
end
push!(jumps, ConstantRateJump(rate11,affect11!))
push!(jumptovars, [6]) # (affect11 changes u[6])


# function rate12(u,p,t)
#     #f_TR ?* DNA_int
#     p[12]*u[6]
# end
# function affect12!(integrator)
#     integrator.u[7] += 1
# end
# push!(jumps, ConstantRateJump(rate12,affect12!))
# push!(jumptovars, [7])


# function rate13(u,p,t)
#     #k_ssRNA_g ?* mRNA_g
#     p[20]*u[7]
# end
# function affect13!(integrator)
#     integrator.u[7] -= 1
#     integrator.u[8] += 1
# end
# push!(jumps, ConstantRateJump(rate13,affect13!))
# push!(jumptovars, [7, 8])

# function rate14(u,p,t)
#     #k_eRNA_g ?* mRNA_g
#     p[17]*u[7]
# end
# function affect14!(integrator)
#     integrator.u[7] -= 1
#     integrator.u[10] += 1
# end
# push!(jumps, ConstantRateJump(rate14,affect14!))
# push!(jumptovars, [7, 10])


# function rate15(u,p,t)
#     #d_RNA_g* mRNA_g
#     p[22]*u[7]
# end
# function affect15!(integrator)
#     integrator.u[7] -= 1
# end
# push!(jumps, ConstantRateJump(rate15,affect15!))
# push!(jumptovars, [7])


# function rate16(u,p,t)
#     #k_dsRNA_ss* mRNA_ss
#     p[21]*u[8]
# end
# function affect16!(integrator)
#     integrator.u[8] -= 1
#     integrator.u[9] += 1
# end
# push!(jumps, ConstantRateJump(rate16,affect16!))
# push!(jumptovars, [8, 9])


# function rate17(u,p,t)
#     #k_eRNA_ss * mRNA_ss
#     p[21]*u[8]
# end
# function affect17!(integrator)
#     integrator.u[8] -= 1
#     integrator.u[11] += 1
# end
# push!(jumps, ConstantRateJump(rate17,affect17!))
# push!(jumptovars, [8, 11])


# function rate18(u,p,t)
#     #d_RNA_ss * mRNA_ss
#     p[23]*u[8]
# end
# function affect18!(integrator)
#     integrator.u[8] -= 1
# end
# push!(jumps, ConstantRateJump(rate18,affect18!))
# push!(jumptovars, [8])


# function rate19(u,p,t)
#     #k_eRNA_ds * mRNA_ds
#     p[19]*u[9]
# end
# function affect19!(integrator)
#     integrator.u[9] -= 1
#     integrator.u[12] += 1
# end
# push!(jumps, ConstantRateJump(rate19,affect19!))
# push!(jumptovars, [9, 12])


# function rate20(u,p,t)
#     #d_RNA_ds * mRNA_ds
#     p[24]*u[9]
# end
# function affect20!(integrator)
#     integrator.u[9] -= 1
# end
# push!(jumps, ConstantRateJump(rate20,affect20!))
# push!(jumptovars, [9])


# function rate21(u,p,t)
#     #k_tpRNA * mRNAc_g
#     p[16]*u[10]
# end
# function affect21!(integrator)
#     integrator.u[10] -= 1
#     integrator.u[18] += 1
# end
# push!(jumps, ConstantRateJump(rate21,affect21!))
# push!(jumptovars, [10, 18])


# function rate22(u,p,t)
#     #d_RNA_g * mRNAc_g
#     p[22]*u[10]
# end
# function affect22!(integrator)
#     integrator.u[10] -= 1
# end
# push!(jumps, ConstantRateJump(rate22,affect22!))
# push!(jumptovars, [10])


# function rate23(u,p,t)
#     #d_RNA_ss * mRNAc_ss
#     p[23]*u[11]
# end
# function affect23!(integrator)
#     integrator.u[11] -= 1
# end
# push!(jumps, ConstantRateJump(rate23,affect23!))
# push!(jumptovars, [11])


# function rate24(u,p,t)
#     #d_RNA_ds * mRNAc_ds
#     p[24]*u[12]
# end
# function affect24!(integrator)
#     integrator.u[12] -= 1
# end
# push!(jumps, ConstantRateJump(rate24,affect24!))
# push!(jumptovars, [12])


# function rate25(u,p,t)
#     #k_trans_f_gGagPol * mRNAc_g
#     p[30]*u[10]
# end
# function affect25!(integrator)
#     integrator.u[13] += 1
# end
# push!(jumps, ConstantRateJump(rate25,affect25!))
# push!(jumptovars, [13])


# function rate26(u,p,t)
#     #k_tpGagPol * P_Tat
#     p[39]*u[13]
# end
# function affect26!(integrator)
#     integrator.u[13] -= 1
#     integrator.u[19] += 1
# end
# push!(jumps, ConstantRateJump(rate26,affect26!))
# push!(jumptovars, [13,19])


# function rate27(u,p,t)
#     #d_pGagPol * P_Tat
#     p[27]*u[13]
# end
# function affect27!(integrator)
#     integrator.u[13] -= 1
# end
# push!(jumps, ConstantRateJump(rate27,affect27!))
# push!(jumptovars, [13])


# function rate28(u,p,t)
#     #k_trans_f_gGag * mRNAc_g
#     p[31]*u[10]
# end
# function affect28!(integrator)
#     integrator.u[14] += 1
# end
# push!(jumps, ConstantRateJump(rate28,affect28!))
# push!(jumptovars, [14])


# function rate29(u,p,t)
#     #k_tpGag * P_Rev
#     p[31]*u[14]
# end
# function affect29!(integrator)
#     integrator.u[14] -= 1
#     integrator.u[20] += 1
# end
# push!(jumps, ConstantRateJump(rate29,affect29!))
# push!(jumptovars, [14, 20])


# function rate30(u,p,t)
#     #d_pGag * P_Rev
#     p[26]*u[14]
# end
# function affect30!(integrator)
#     integrator.u[14] -= 1
# end
# push!(jumps, ConstantRateJump(rate30,affect30!))
# push!(jumptovars, [14])



# function rate31(u,p,t)
#     #k_trans_f_ssgp160 * mRNAc_ss
#     p[32]*u[11]
# end
# function affect31!(integrator)
#     integrator.u[15] += 1
# end
# push!(jumps, ConstantRateJump(rate31,affect31!))
# push!(jumptovars, [15])


# function rate32(u,p,t)
#     #k_trans_f_ssgp160 * P_GagPol
#     p[32]*u[15]
# end
# function affect32!(integrator)
#     integrator.u[15] -= 1
#     integrator.u[21] += 1
# end
# push!(jumps, ConstantRateJump(rate32,affect32!))
# push!(jumptovars, [15, 21])


# function rate33(u,p,t)
#     #d_pgp160 * P_GagPol
#     p[25]*u[15]
# end
# function affect33!(integrator)
#     integrator.u[15] -= 1
# end
# push!(jumps, ConstantRateJump(rate33,affect33!))
# push!(jumptovars, [15])


# function rate34(u,p,t)
#     #k_trans_f_dsTat
#     p[33]*u[12]
# end
# function affect34!(integrator)
#     integrator.u[16] += 1
# end
# push!(jumps, ConstantRateJump(rate34,affect34!))
# push!(jumptovars, [16])


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

function inhibition(D, IC50, n)
    return 1 / ((D/IC50)^n + 1)
end

function apply_inhibitors!(p, inhibitors, conc2, conc3)
    concentrations = [conc2, conc3]

    eff_rt = 1.0
    # eff_mat = 1.0

    for (i, row) in enumerate(eachrow(inhibitors))
        eff = inhibition(concentrations[i], row.IC50, row.m)
        # if row.target == "k_mat"
        #     eff_mat *= eff
        if  row.target == "k_RT"
            eff_rt *= eff
        end
    end
    p[4] *= eff_rt
    # p[49] *= eff_mat
end

# parameters
pt = CSV.File("parameters.txt", header=false, delim=' ', ignorerepeated=true)
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
matrix = zeros(n_NRTI, n_NRTI)

for (idx2, conc2) in enumerate(NRTI_conc_range)
    for (idx3, conc3) in enumerate(NRTI_conc_range)
                
        local_p = copy(p)
        apply_inhibitors!(local_p, inhibitors, conc2, conc3) 
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
        matrix[idx2, idx3] = local_p_infected / coeff
    end
end

println(size(matrix))

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

h5open("data/matrix_u1.h5", "w") do file
    write(file, "matrix", matrix)
end

println("Done!")