#' Counts the frequency of different aa
#'
#' @param path_in input fasta file
#'
#' @return data frame
get_aa_composition = function(path_in){

    x = CHNOSZ::read.fasta(path_in)
    l = names(ape::read.FASTA(path_in))
    n = base::strsplit(l, "|", fixed = TRUE)
    s = data.frame("entry_code" = sapply(n, `[`, 1), "organism" = sapply(n, `[`, 2))

    x = x[ , ! colnames(x) %in% c("protein", "organism", "ref", "abbrv", "chains") ]

    data.frame(s, x)
}


#' Counts the frequency of chemical elements in a sequence
#'
#' Duplicated code from get_aa_composition()!!! Element counting depends on the
#' aa count from CHNOSZ::read.fasta, but get_aa_composition saves a data.frame
#'
#'
#' @param path_in input fasta file
#'
#' @return data frame
get_element_composition = function(path_in){

    x = CHNOSZ::read.fasta(path_in)
    y = CHNOSZ::protein.formula(x)
    l = names(ape::read.FASTA(path_in))
    n = base::strsplit(l, "|", fixed = TRUE)
    s = data.frame("entry_code" = sapply(n, `[`, 1),
                   "organism"   = sapply(n, `[`, 2))

    data.frame(s, y, row.names = NULL)
}

