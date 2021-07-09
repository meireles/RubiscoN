library("targets")
library("corrplot")
library("ppcor")

# load data
tar_load(final_data)

# Run correlations for each clade in final_data

log_try                 = TRUE
normalize_by_seq_length = TRUE

for(i in names(final_data)){

    message(i)

    try = final_data[[i]][["try_data"]]
    aa  = final_data[[i]][["aa_comp"]]
    ele = final_data[[i]][["elem_comp"]]

    # pick columns of interest
    col_keep_try = c("LfArea_Lflt", "SLA", "N_Dry", "P_Dry") # "VcM_Dry"
    col_keep_aa  = ! colnames(aa) %in% c("entry_code", "organism" )
    col_keep_ele = ! colnames(ele) %in% c("entry_code", "organism" )

    # subset datasets
    try = try[ , col_keep_try]
    aa  = aa[ , col_keep_aa]
    ele = ele[ , col_keep_ele]


    # Log TRY data
    if(log_try){
        try = log(try)
    }

    # Normalize by sequence length?
    if(normalize_by_seq_length){
        s   = rowSums(aa)
        aa  = aa / s
        ele = ele / s
    }

    # add N : C to the element matrix
    ele = cbind(ele, "N:C" = ele$N / ele$C)

    # combine data for correlations
    mat_aa_try = cbind(aa , try)
    mat_ele_try = cbind(ele , try)

    mat_aa_try  = na.omit(mat_aa_try)
    mat_ele_try = na.omit(mat_ele_try)


    message(nrow(mat_ele_try))

    ####################
    # regular correlation
    ####################

    # corr traits / aa freq
    cmat_aa       = cor(mat_aa_try)
    cmat_aa_p     = cor.mtest(mat_aa_try)$p

    # corr traits / element freq
    cmat_ele      = cor(mat_ele_try)
    cmat_ele_p    = cor.mtest(mat_ele_try)$p

    ####################
    # partial correlation
    ####################

    p_c_aa          = pcor(mat_aa_try)
    p_cmat_aa       = p_c_aa$estimate
    p_cmat_aa_p     = p_c_aa$p.value

    dimnames(p_cmat_aa)   = list(colnames(mat_aa_try), colnames(mat_aa_try))
    dimnames(p_cmat_aa_p) = list(colnames(mat_aa_try), colnames(mat_aa_try))

    # corr traits / element freq
    p_c_ele         = pcor(mat_ele_try)
    p_cmat_ele      = p_c_ele$estimate
    p_cmat_ele_p    = p_c_ele$p.value

    dimnames(p_cmat_ele)   = list(colnames(mat_ele_try), colnames(mat_ele_try))
    dimnames(p_cmat_ele_p) = list(colnames(mat_ele_try), colnames(mat_ele_try))

    # Plot
    pdf(file = paste0("kitchensink/tmp_figs/", "cor_try_rubisco_", i, ".pdf"),
        width = 10, height = 12)

    par(mfrow = c(2, 2), oma = c(2, 2, 10, 2))


    corrplot(corr = cmat_aa, type = "upper", p.mat = cmat_aa_p,
             sig.level = 0.005, insig = "blank", diag = FALSE,
             addCoef.col = "black", addCoefasPercent = TRUE, number.cex = 0.4)

    title("Correlation")


    corrplot(corr = p_cmat_aa, type = "upper", p.mat = p_cmat_aa_p,
             sig.level = 0.005, insig = "blank", diag = FALSE,
             addCoef.col = "black", addCoefasPercent = TRUE, number.cex = 0.4)
    title("Partial Correlation")

    corrplot(corr = cmat_ele, type = "upper", p.mat = cmat_ele_p,
             sig.level = 0.005, insig = "blank", diag = FALSE,
             addCoef.col = "black", addCoefasPercent = TRUE, number.cex = 0.4)


    corrplot(corr = p_cmat_ele, type = "upper", p.mat = p_cmat_ele_p,
             sig.level = 0.005, insig = "blank", diag = FALSE,
             addCoef.col = "black", addCoefasPercent = TRUE, number.cex = 0.4)

    mtext(paste0(i, " (N = ", nrow(mat_ele_try), ")"), outer = TRUE, cex = 1.5, line = 3)

    dev.off()
}
