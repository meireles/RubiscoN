library("targets")
library("corrplot")

# load data
tar_load(final_data)


# Run correlations for each clade in final_data

for(i in names(final_data)){

    try = final_data[[i]][["try_data"]]
    aa  = final_data[[i]][["aa_comp"]]
    ele = final_data[[i]][["elem_comp"]]

    # pick columns of interest
    col_keep_try = c("LfArea_Lflt", "SLA", "N_Dry", "P_Dry", "VcM_Dry")
    col_keep_aa  = ! colnames(aa) %in% c("entry_code", "organism" )
    col_keep_ele = ! colnames(ele) %in% c("entry_code", "organism" )

    # subset datasets
    try = try[ , col_keep_try]
    aa  = aa[ , col_keep_aa]
    ele = ele[ , col_keep_ele]

    # combine data for correlations
    mat_aa_try = cbind(aa , try)
    mat_ele_try = cbind(ele , try)


    # Pick use param for cor
    use = "pairwise.complete.obs" # "complete.obs"


    # corr traits / aa freq
    cmat_aa       = cor(mat_aa_try, use = use)
    cmat_aa_p     = cor.mtest(mat_aa_try, use = use)

    # corr traits / element freq
    cmat_ele      = cor(mat_ele_try, use = use)
    cmat_ele_p    = cor.mtest(mat_ele_try, use = use)


    # Plot
    pdf(file = paste0("kitchensink/tmp_figs/", "cor_try_rubisco_", i, ".pdf"),
        width = 6, height = 12)

    par(mfrow = c(2, 1))

    corrplot(corr = cmat_aa, type = "upper", p.mat = cmat_aa_p$p,
             sig.level = 0.005, insig = "blank", diag = FALSE,
             addCoef.col = "black", addCoefasPercent = TRUE, number.cex = 0.4)

    corrplot(corr = cmat_ele, type = "upper", p.mat = cmat_ele_p$p,
             sig.level = 0.005, insig = "blank", diag = FALSE,
             addCoef.col = "black", addCoefasPercent = TRUE, number.cex = 0.4)

    mtext(i, outer = TRUE, cex = 1.5)

    dev.off()
}

