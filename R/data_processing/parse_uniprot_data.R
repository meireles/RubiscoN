#' Parse sequence from the raw uniprot data
#'
#' @param x data downloaded using `download_uniprot_raw_data`
#' @param path_out save fasta file to
#'
#' @return returns path_out (required by targets) and side effect of saving data
parse_uniprot_fasta = function(x, path_out){
    cat(paste0(">", x$Entry,"|", x$Organism,"\n", x$Sequence,"\n", collapse = ""),
        file = path_out, append = FALSE)
    path_out
}


#' Parse uniprot features
#'
#' @param x data downloaded using `download_uniprot_raw_data`
#' @param path_out save fasta file to
#' @param feature Features to retrieve
#'
#' @return returns path_out (required by targets) and side effect of saving data
parse_uniprot_features = function(x,
                                  path_out,
                                  feature = c("Active site"   = "ACT_SITE",
                                              "Binding site"  = "BINDING",
                                              "Metal binding" = "METAL",
                                              "Site"          = "SITE")){

    r = setNames(vector("list", nrow(x)), x$Entry)

    pb    = txtProgressBar(min = 0, max = length(r), initial = 0, style = 3)

    for(i in seq_len(nrow(x))){

        y = x[i, ]

        s = setNames(vector("list", length(feature)), names(feature))

        for(j in seq_along(feature)){
            a = y[[names(feature[j])]]
            b = base::strsplit(a, paste0(feature[j], " |/note=|/evidence="))
            c = lapply(b, function(z){
                z1 = z[z != ""]
                z2 = gsub(";|\"", "", z1)
                z3 = str_trim(z2)
                z4 = matrix(z3, ncol = 3, byrow = TRUE,
                            dimnames = list(NULL, c("position", paste0(names(feature[j]), " feature"), "")))
                as.data.frame(z4[ , 1:2, drop = FALSE])
            })
            s[[names(feature[j])]] = c[[1]]
        }

        r[[i]] = Reduce(function(x1, x2){
            merge(x1, x2, by = "position", all = TRUE)}, s)

        r[[i]] = data.frame("Entry" = y[ , "Entry"], r[[i]], check.names = FALSE)

        setTxtProgressBar(pb = pb, value = i)
    }
    close(pb)

    s = do.call(rbind, r)
    write_tsv(s, path_out)
    path_out
}
