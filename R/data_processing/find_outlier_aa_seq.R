#' Find aa sequences that are too different from the rest
#'
#' @param z sequences of class DNAbin (from ape)
#' @param k_vec vector of kmers to compute distances with
#' @param sigma sequences more distant than `sigma` standard deviations are
#'              considered too different and excluded
#'
#' @return DNAbin sequences
find_outlier_aa_seq = function(z, k_vec = c(4, 5), sigma = 6){
    m = matrix(data     = 0,
               nrow     = length(z),
               ncol     = length(k_vec) + 1,
               dimnames = list(names(z), c(k_vec, "outlier_n_hits") ))

    for(i in seq_along(k_vec)){
        message("kmer k = ", k_vec[i])
        y = kmer::kcount(z, k = k_vec[i], residues = "AA")
        x = t( t(y) - colMeans(y) )
        q = rowSums(abs(x))
        p = log(q)

        m[ , i] = p

        l = mean(p) - sigma * sd(p)
        h = mean(p) + sigma * sd(p)
        o = p > h | p < l

        m[ , "outlier_n_hits"] = m[ , "outlier_n_hits"] + o
    }
    m
}
