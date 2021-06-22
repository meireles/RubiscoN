#' Align with MAFFT
#'
#' @param path_in path for the unaligned data
#' @param path_out save fasta alignment to
#'
#' @return returns path_out (needed by targets). Saves fasta as a side effect
mafft = function(path_in, path_out){
    command = paste("bin/mafft-mac/mafft.bat",
                    "--thread 8 --threadtb 2 --threadit 0",
                    "--reorder --large --retree 4 --op 2.34 --ep 0.123",
                    path_in, ">", path_out, sep = " ")

    base::system(command)

    path_out
}
