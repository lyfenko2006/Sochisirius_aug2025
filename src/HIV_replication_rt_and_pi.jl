using DifferentialEquations
using PyPlot
using LaTeXStrings
using Printf
using CSV
using JSON
using HDF5
using DataFrames
using DelimitedFiles
using JLD2

ENV["JULIA_SSL_LIBRARY"] = "/usr/lib/x86_64-linux-gnu/libssl.so.3"  
using Pkg
Pkg.build("LibCURL")  # Пересобираем LibCURL

function HIV_replication(du,u,p,t)
    V_free, V_bound, RNA_cor, DNA_cor, DNA_nuc, DNA_int,
    	mRNA_g, mRNA_ss, mRNA_ds, mRNAc_g, mRNAc_ss, mRNAc_ds,
    	P_Tat, P_Rev, P_GagPol, P_Gag, P_gp160,
    	P_memGagPol, P_memGag, P_memgp160, RNA_mem,
    	V_previrion, V_bud, V_mat = u

    #@unpack_Params p
    (k_bound, d_bound, k_fuse,   	  	        # binding & fusion
     k_RT, d_RNA_cor, d_DNA_cor, 	  	        # reverse transcription
     k_DNA_t, d_DNA_nuc, k_int, d_DNA_int, 	        # integration
     TR_cell, TR_Tat, theta_Tat, theta_Rev, beta, 	# transcription
     k_tpRNA, k_eRNA_g, k_eRNA_ss, k_eRNA_ds, 	        # export and transport
     k_ssRNA_g, k_dsRNA_ss,			        # splicing
     d_RNA_g, d_RNA_ss, d_RNA_ds,		        # \
     d_pgp160, d_pGag, d_pGagPol,       	        #  turnover
     d_pTat, d_pRev, 				        # /
     k_trans_f_gGagPol, k_trans_f_gGag,		        # translation rates of ...
     k_trans_f_ssgp160, k_trans_f_dsTat,		# the proteins Gag-Pol, Gag, gp160, ...
     k_trans_f_dsRev,					# Tat, Rev
     d_memGagPol, d_memGag, d_memgp160,                 # degradation at membrane
     k_tpGag, k_tpGagPol, k_tpgp160,                    # transport to membrane
     k_comb, k_comb_N_RNA, k_comb_N_Gag,                # assembly
     k_comb_N_GagPol, k_comb_N_gp160, d_comb,	        # //----//
     k_bud, d_bud, k_mat, d) = p		        # budding & maturation


    #dV_free
    du[1] = -k_bound*V_free - d*V_free

    #dV_bound
    du[2] = k_bound*V_free - (k_fuse+d_bound)*V_bound

    #dRNA_cor
    du[3] = k_fuse*V_bound - (k_RT+d_RNA_cor)*RNA_cor

    #dDNA_cor
    du[4] = k_RT*RNA_cor - (k_DNA_t+d_DNA_cor)*DNA_cor

    #dDNA_nuc
    du[5] = k_DNA_t*DNA_cor - (k_int+d_DNA_nuc)*DNA_nuc

    #dDNA_int
    du[6] = k_int*DNA_nuc - d_DNA_int*DNA_int

    #f_Tat, f_Rev, TR
    f_Tat = P_Tat/(P_Tat + theta_Tat)
    f_Rev = P_Rev/(P_Rev + theta_Rev)
    TR = TR_cell + f_Tat*TR_Tat

    #dmRNA_g
    du[7] = TR*DNA_int - (k_eRNA_g*f_Rev + k_ssRNA_g*(1.0-beta*f_Rev) + d_RNA_g)*mRNA_g

    #dmRNA_ss
    du[8] = (1.0-beta*f_Rev)*k_ssRNA_g*mRNA_g - (k_eRNA_ss*f_Rev+d_RNA_ss+k_dsRNA_ss*(1.0-beta*f_Rev))*mRNA_ss

    #dmRNA_ds
    du[9] = (1.0-beta*f_Rev)*k_dsRNA_ss*mRNA_ss - (d_RNA_ds+k_eRNA_ds)*mRNA_ds

    #dmRNAc_g
    du[10] = f_Rev*k_eRNA_g*mRNA_g - (k_tpRNA + d_RNA_g)*mRNAc_g

    #dmRNAc_ss
    du[11] = f_Rev*k_eRNA_ss*mRNA_ss - d_RNA_ss*mRNAc_ss

    #dmRNAc_ds
    du[12] = k_eRNA_ds*mRNA_ds - d_RNA_ds*mRNAc_ds

    #dP_Tat
    du[13] = k_trans_f_dsTat*mRNAc_ds - d_pTat*P_Tat

    #dP_Rev
    du[14] = k_trans_f_dsRev*mRNAc_ds - d_pRev*P_Rev

    #dP_GagPol
    du[15] = k_trans_f_gGagPol*mRNAc_g - (k_tpGagPol+d_pGagPol)*P_GagPol

    #dP_Gag
    du[16] = k_trans_f_gGag*mRNAc_g - (k_tpGag+d_pGag)*P_Gag

    #P_gp160
    du[17] = k_trans_f_ssgp160*mRNAc_ss - (k_tpgp160+d_pgp160)*P_gp160

    #dP_memGagPol
    du[18] = k_tpGagPol*P_GagPol - k_comb_N_GagPol*RNA_mem*P_memGag*P_memGagPol*P_memgp160 - d_memGagPol*P_memGagPol

    #dP_memGag
    du[19] = k_tpGag*P_Gag - k_comb_N_Gag*RNA_mem*P_memGag*P_memGagPol*P_memgp160 - d_memGag*P_memGag

    #dP_memgp160
    du[20] = k_tpgp160*P_gp160 - k_comb_N_gp160*RNA_mem*P_memGag*P_memGagPol*P_memgp160 - d_memgp160*P_memgp160

    #dRNA_mem
    du[21] = k_tpRNA*mRNAc_g - k_comb_N_RNA*RNA_mem*P_memGag*P_memGagPol*P_memgp160 - d_RNA_g*RNA_mem

    #dV_previrion
    du[22] = k_comb*RNA_mem*P_memGag*P_memGagPol*P_memgp160 - (k_bud+d_comb)*V_previrion

    #dV_bud
    du[23] = k_bud*V_previrion - (k_mat+d_bud)*V_bud

    #dV_mat
    du[24] = k_mat*V_bud - d*V_mat

    #cumulative viral load
    du[25] = V_mat

    du[26] = u[12]

    du[27] = k_mat * u[23]
end

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

function inhibition(D, IC50, n)
    return 1 / ((D/IC50)^n + 1)
end

function apply_inhibitors!(p, inhibitors, conc1, conc2, conc3)
    concentrations = [conc1, conc2, conc3]

    eff_rt = 1.0
    eff_mat = 1.0

    for (i, row) in enumerate(eachrow(inhibitors))
        eff = inhibition(concentrations[i], row.IC50, row.m)
        if row.target == "k_mat"
            eff_mat *= eff
        elseif  row.target == "k_RT"
            eff_rt *= eff
        end
    end
    p[4] *= eff_rt
    p[49] *= eff_mat
end

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
            apply_inhibitors!(local_p, inhibitors, conc1, conc2, conc3)
            
            local_prob = ODEProblem(HIV_replication,local_u0,tspan,local_p)
            local_sol = solve(local_prob,RadauIIA5(),reltol=1e-14,abstol=1e-14)

            local_u27 = [u[27] for u in local_sol.u]
            matrix[i, j, k] = local_u27[length(local_u27)] / coeff
        end
    end
end

# h5open("data/matrix_u2.h5", "w") do file
#     write(file, "matrix", matrix)
# end
data_dict = Dict("matrix" => matrix)
safesave(datadir("sims", "matrix_u2.jld2"), data_dict)
println("Done!")

# # draw numerical solution
# L = [L"V_{free}", L"V_{bound}", L"RNA_{cor}", L"DNA_{cor}", L"DNA_{nuc}", L"DNA_{int}",
#     L"mRNA_g", L"mRNA_{ss}", L"mRNA_{ds}", L"mRNAc_g", L"mRNAc_{ss}", L"mRNAc_{ds}",
#     L"P_{Tat}", L"P_{Rev}", L"P_{Gag-Pol}", L"P_{Gag}", L"P_{gp160}",
#     L"P_{mem,Gag-Pol}", L"P_{mem,Gag}", L"P_{mem,gp160}", L"RNA_{mem}",
#     L"V_{pre-virion}", L"V_{bud}", L"V_{mat}", L"\int_0^{t} V_{mat} dt"]
# i2mm = 0.039370077777777776



# function draw_sol(sol, filename; ts=nothing, p24_ts=nothing)
#     function mydraw(disp_vars, tspan; loc="upper left", ythresh=1.0, ythscale=0.2,
#                                       lstyle="-", colors=nothing)
#         if isnothing(colors)
#             pl = ax.plot(sol.t, Array(sol')[:,disp_vars], linestyle=lstyle)
#         else
#             for (var,col) in zip(disp_vars,colors)
#                 pl = ax.plot(sol.t, Array(sol')[:,var], linestyle=lstyle, color=col)
#             end
#         end
#         ax.set_xlim(tspan)
#         ax.set_xlabel("t, hours")
#         #ax.set_ylim(bottom=0)
#         ax.set_ylim((0,1.15*ax.get_ylim()[2]))
#         ax.set_yscale("symlog", linthresh=ythresh, linscale=ythscale)
#         #ax.tick_params(axis="y", which="minor", bottom=false)
#         #tick_params(axis="y", which="minor", left=true)
#         #ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
#         ax.grid(which="major", linestyle="-")
#         ax.minorticks_on()
#         ax.legend(loc=loc, labels=L[disp_vars])
#         #y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=false)
#         #y_formatter.set_powerlimits((-1,0.5))
#         #y_formatter.set_powerlimits((-1,100))
#         #ax.yaxis.set_major_formatter(y_formatter)
#         return pl
#     end

#     PyPlot.ioff()
#     try
#         PyPlot.rc("text", usetex=true)
#     catch e
#         @warn "Latex not available, falling back to non-Latex text"
#         PyPlot.rc("text", usetex=false)
#     end

#     PyPlot.rc("font", size=12)

#     fig, axes = subplots(3,3,
#          figsize=(270.0*i2mm,210.0*i2mm),
#          constrained_layout=true)

#     ax = axes[1,1]
#     mydraw([1,2], [0,5], loc="upper right", ythresh=1e-2)
#     ax.set_yscale("linear")
#     ax.set_ylabel("particles")
#     ax.set_title("Binding and fusion")

#     ax = axes[1,2]
#     mydraw([3,4,5,6], [0,36], loc="center right", ythresh=1e-4)
#     ax.set_yscale("linear")
#     ax.set_ylabel("molecules")
#     ax.set_title("RT and integration")

#     ax = axes[1,3]
#     pl = mydraw([7,8,9,], [0,36], lstyle="--")
#     mydraw([10,11,12], [0,36], colors=[l.get_color() for l in pl])
#     ax.set_ylim((0,1.15*maximum(Array(sol')[:,[7,8,9,10,11,12]])))
#     ax.legend(loc="lower right", labels=L[[7,8,9,10,11,12]])
#     ax.set_ylabel("molecules")
#     ax.set_title("Transcription")

#     ax = axes[2,1]
#     mydraw([13,14], [0,36], loc="lower right")
#     ax.set_ylabel("molecules")
#     ax.set_title("Regulatory proteins")

#     ax = axes[2,2]
#     mydraw([15,16,17], [0,36], loc="lower right")
#     ax.set_ylabel("molecules")
#     ax.set_title("Protein translation")

#     ax = axes[2,3]
#     mydraw([18,19,20,21], [0,36])
#     ax.set_ylabel("molecules")
#     ax.set_title("Proteins at membrane")

#     ax = axes[3,1]
#     mydraw([22,23], [0,36])
#     ax.set_ylabel("particles")
#     ax.set_title("Assembly and budding")

#     ax = axes[3,2]
#     mydraw([24], [0,36])
#     ax.set_yscale("linear")
#     ax.set_ylabel("particles")
#     ax.set_title("Maturation")

#     if !isnothing(ts) && !isnothing(p24_ts)
#         ax = axes[3,3]
#         ax.plot(sol.t, Array(sol')[:,24])
#         ax.scatter(ts, p24_ts, marker="+", color="k")
#         ax.set_xlim([0,36])
#         ax.set_xlabel("t, hours")
#         ax.set_ylim((0,1.15*ax.get_ylim()[2]))
#         ax.grid(which="major", linestyle="-")
#         ax.minorticks_on()
#         ax.legend(loc="upper left", labels=[L[24], "exp. data"])
#         ax.set_yscale("linear")
#         ax.set_ylabel("particles")
#         ax.set_title("Viral release")
#     else
#         ax = axes[3,3]
#         mydraw([25], [0,36])
#         ax.set_yscale("linear")
#         ax.set_ylabel("particles")
#         ax.set_title("Cumulative number of released virions")
#     end

#     try
#         savefig(filename)
#     catch e 
#         @warn "Failed to save as $filename, trying PNG"
#         newname = replace(filename, r"\.pdf$" => ".png")
#         savefig(newname)
#     end

#     plt.close(fig)
# end
# draw_sol(sol, "plot_model_3TC_AZT.pdf")
