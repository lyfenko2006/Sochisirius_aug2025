module All_function

using Revise
using DrWatson
using DifferentialEquations
using PyPlot
using LaTeXStrings
using Printf
using CSV
using JSON
using HDF5
using DataFrames
using DelimitedFiles

ENV["JULIA_SSL_LIBRARY"] = "/usr/lib/x86_64-linux-gnu/libssl.so.3.3"  
using Pkg
Pkg.build("LibCURL")  # Пересобираем LibCURL
export HIV_replication, inhibition, apply_inhibitors_1!, apply_inhibitors_2!
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

function inhibition(D, IC50, n)
    return 1 / ((D/IC50)^n + 1)
end

function apply_inhibitors_1!(p, inhibitors, conc1)
    concentrations = [conc1]

    eff_mat = 1.0

    for (i, row) in enumerate(eachrow(inhibitors))
        eff = inhibition(concentrations[i], row.IC50, row.m)
        if row.target == "k_mat"
            eff_mat *= eff
        end
    end

    p[49] *= eff_mat
end

function apply_inhibitors_2!(p, inhibitors, conc1, conc2, conc3)
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

end