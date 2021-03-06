#' Fist clean up steps for the aa data.
#'
#' removes outlier, uncertain, and short seqs
#'
#' @param path_in parth to imput fasta file
#' @param path_out save fasta file to
#' @param uncertain_tol percentage of uncertain aa deemed as too many
#' @param k_vec vector of kmers to compute distances with
#' @param sigma sequences more distant than `sigma` standard deviations are
#'              considered too different and excluded
#' @param min_length only keep seqs longer than `min_length`
#'
#' @return returns path_out (needed by targets). Saves fasta as a side effect
initial_seq_cleanup = function(path_in,
                               path_out,
                               uncertain_tol = 0.05,
                               k_vec         = c(4, 5),
                               sigma         = 6,
                               min_length){

    x = ape::read.FASTA(path_in, "AA")

    message("Finding uncertain sequences...")
    y = remove_uncertain_seqs(x, uncertain_tol = uncertain_tol)


    message("Removing sequences shorter than 400...")
    k = sapply(y, length) >= min_length
    y = y[k]

    message("Finding outlier sequences...")
    o = find_outlier_aa_seq(y, sigma = sigma, k_vec = k_vec)

    z = y[ o[ , "outlier_n_hits"] < length(k_vec) ]

    ape::write.FASTA(z, file = path_out)

    path_out
}


trim_alignment_ends = function(path_in, path_out, prop_gap = 0.1){
    x = Biostrings::readAAStringSet(path_in, format = "fasta")
    y = as.matrix(x)
    z = apply(y, 2, function(x){ sum(x == "-") })
    r = range(which(z <= max(z) * prop_gap))
    k = seq.int(r[1], r[2] )
    a = ape::as.AAbin( y[ , k])

    ape::write.FASTA(a, path_out, append = FALSE)
    path_out
}


clean_tax_and_pick_long_seq = function(path_in, path_out){

    x = Biostrings::readAAStringSet(path_in, format = "fasta")
    y = as.matrix(x)
    n = sapply(base::strsplit(rownames(y), "|", fixed = TRUE), `[`, 2)
    g = which(n == "NA")
    y = y[ -g, ]
    n = n[ -g ]
    d = sort(unique(n[duplicated(n)]))

    message("Keeping longest sequence from duplicated species...")

    r = NULL

    for(i in d){
        message(i)
        w = which(n %in% i)
        q = apply( y[w, ], 1, function(h){ sum(h == "-") })

        r = c(r, w[ - which.max(q) ]) #yes, I am ashamed of this line
    }

    y = y[ - r, ]
    n = n[ - r ]
    y = y[ - which_genera_and_duplicated_species(n), ]

    ape::write.FASTA(ape::as.AAbin(y), path_out, append = FALSE)
    path_out
}

