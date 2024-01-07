#!/usr/bin/env

##Check if all the positional variables need are provided, if not then exits and prints out the requirements

.libPaths('/lustre/scratch126/tol/teams/lawniczak/users/jr35/local_libs/R410_libs')

##load libraries
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(BiocParallel)
library(purrr)

##Initializes the positional variables
args = commandArgs(trailingOnly=TRUE)

i <- args[1]

msc_gcombi_noLE_sub_sce <- readRDS('/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/pseudotime_tradeseq_inp/msc_gcombi_noLE_sub_sce.RDS')
# msc_gcombi_noLE_sub_sds <- readRDS('/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/pseudotime_tradeseq_inp/msc_gcombi_noLE_sub_sds.RDS')
msc_gcombi_noLE_sub_ptime <- readRDS('/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/pseudotime_tradeseq_inp/msc_gcombi_noLE_sub_ptime.RDS')
msc_gcombi_noLE_sub_cw <- readRDS('/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/pseudotime_tradeseq_inp/msc_gcombi_noLE_sub_cw.RDS')

# condtns <- c('strain_dset','strain_dset','strain_dset','strain_dset','strain_dset')
# condtns <- c('dset','dset','dset','dset','dset','strain_dset','strain_dset','strain_dset','strain_dset','strain_dset' ,'strain','strain','strain','strain','strain')
condtns <- c(rep(c('dset'), 8), rep(c('lf_strain'), 20))
# comps <- c("msc14_f", "msc14_m", "msc1_f","msc1_m", "msc3_f", "msc3_m", "msc13_f", "msc13_m", "msc14_f", "msc14_m")
comps <- names(msc_gcombi_noLE_sub_sce)
print(comps)

##fitting GAMs
knots_no <- c(6,4,6,4,6,4,4,4,4,4,4,4,4,4,4,4,4,4,6,6,4,4,4,4,6,6,4,4)
# knots_no <- c(8,7,8,9,8,6,6,6,8,8, 6,6,6,6,6)

# sce_obj = c(msc_gcombi_noLE_sub_sce, msc_gcombi_noLE_sub_sce, msc_gcombi_noLE_sub_sce)
# sce_obj = c(msc_gcombi_noLE_sub_sce,msc_gcombi_noLE_sub_sce)
sce_obj = msc_gcombi_noLE_sub_sce
# sce_obj = c(msc_gcombi_noLE_sub_sce[1:5], msc_gcombi_noLE_sub_sce[1:5],msc_gcombi_noLE_sub_sce[6:10], msc_gcombi_noLE_sub_sce[6:10],msc_gcombi_noLE_sub_sce[c(11:15)])

# ptime_obj = c(msc_gcombi_noLE_sub_ptime, msc_gcombi_noLE_sub_ptime, msc_gcombi_noLE_sub_ptime)
# ptime_obj = c(msc_gcombi_noLE_sub_ptime, msc_gcombi_noLE_sub_ptime)
ptime_obj = msc_gcombi_noLE_sub_ptime
# ptime_obj = c(msc_gcombi_noLE_sub_ptime[1:5], msc_gcombi_noLE_sub_ptime[1:5], msc_gcombi_noLE_sub_ptime[6:10], msc_gcombi_noLE_sub_ptime[6:10], msc_gcombi_noLE_sub_ptime[c(11:15)])
# cw_obj = c(msc_gcombi_noLE_sub_cw, msc_gcombi_noLE_sub_cw, msc_gcombi_noLE_sub_cw)
# cw_obj = c(msc_gcombi_noLE_sub_cw, msc_gcombi_noLE_sub_cw)

cw_obj = msc_gcombi_noLE_sub_cw
# cw_obj = c(msc_gcombi_noLE_sub_cw[1:5], msc_gcombi_noLE_sub_cw[1:5], msc_gcombi_noLE_sub_cw[6:10], msc_gcombi_noLE_sub_cw[6:10], msc_gcombi_noLE_sub_cw[c(11:15)])

i=as.numeric(i)

## Subset the relevant cells for each run
substd_sce = sce_obj[[i]][,!is.na(sce_obj[[i]]@colData[,condtns[i]])]
substd_ptime = ptime_obj[[i]][!is.na(sce_obj[[i]]@colData[,condtns[i]])]
substd_cw = cw_obj[[i]][!is.na(sce_obj[[i]]@colData[,condtns[i]])]

# print(str(substd_sce))
# print(condtns[i])
# print(comps[i])

## Tradeseq script
msc14_ref_ts_sce_ub_de  <- fitGAM(counts = counts(substd_sce), 
                    pseudotime = substd_ptime, 
                    cellWeights = substd_cw, 
                    nknots =knots_no[i],
                    genes = rownames(substd_sce)[rowSums(counts(substd_sce) != 0) > 5],
                    conditions = factor(SummarizedExperiment::colData(substd_sce)[, condtns[i]]), 
                    verbose = T, 
                    parallel = T
                    )
            
##Bal status
# bst <- c(rep('_unbal', 10), rep('', 10), rep('', 5))
bst <- c(rep('', 28))

saveRDS(msc14_ref_ts_sce_ub_de, paste0('/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/pseudotime_tradeseq_op/tseq_obj_', condtns[i], '_', comps[i],bst[i], '.RDS'))

##Running condition test
msc14_condRes_ub_de <- conditionTest(msc14_ref_ts_sce_ub_de, pairwise = T)
            
saveRDS(msc14_condRes_ub_de, paste0('/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali/data/processed/DB52e_star/pseudotime_tradeseq_op/tseq_cond_obj_', condtns[i], '_', comps[i], bst[i], '.RDS'))
