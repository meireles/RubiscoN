#' Remove sequences with too many uncertain aa
#'
#' @param x sequence of DNAbin class (from ape)
#' @param uncertain_tol percentage of uncertain aa deemed as too many.
#'
#' @return DNAbin seq
remove_uncertain_seqs = function(x, uncertain_tol){
    amb_c = c("X", "B", "Z", "J", "O", "U", "*", ".")
    amb_r = sapply(amb_c, charToRaw)
    p_amb = t(sapply(x, function(z){ sum(z %in% amb_r) / length(z) }))
    keep  = p_amb <= uncertain_tol
    x[keep]
}
