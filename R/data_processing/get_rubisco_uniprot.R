#' Download uniprot data
#'
#' @param path_out save file to
#' @param min_length min aa seq length
#' @param max_length max aa seq length
#' @param gene Name of gene that encodes the protein. Defaults to rbcl
#' @param taxon uniprot taxon id. Defaults to 58023 ("Tracheophyta")
#' @param request_lim max number of entries to request at once.
#'
#' @return returns path_out (required by targets) and side effect of saving data
#' to path_out
download_uniprot_raw_data = function(path_out,
                                     min_length  = 400,
                                     max_length  = 505,
                                     gene        = "rbcl",
                                     taxon       = 58023,
                                     request_lim = 10000){
    ## build request
    request_base = "https://www.uniprot.org/uniprot/?query="
    query        = paste0("taxonomy:", taxon, "+AND+",
                          "gene_exact:", gene, "+AND+",
                          "length:[",
                          min_length, " TO ", max_length,"]")
    format       = "&format=tab"
    base_cols    = "&columns=id,reviewed,organism,length,annotation_score,sequence"
    struct_cols  = ",feature(BETA_STRAND),feature(HELIX),feature(TURN)"
    feature_cols = ",feature(ACTIVE_SITE),feature(BINDING_SITE),feature(DNA_BINDING),feature(METAL_BINDING),feature(NP_BIND),feature(SITE)"

    n_entries_q     = gsub(" ", "%20", paste0(request_base, query, "&format=list"))
    n_entries       = str_count(httr::content(httr::POST(n_entries_q), "text"),  "\n")

    cat(n_entries, " sequences found.")

    request_offset  = request_lim * seq(0, floor(n_entries / request_lim) )

    request_pars    = paste0("&limit=", request_lim,"&offset=", request_offset)
    request         = gsub(" ", "%20", paste0(request_base, query, format,
                                              base_cols, struct_cols, feature_cols,
                                              request_pars))

    # Iterate over number of requests
    for(i in seq_along(request)){
        response   = httr::POST(request[[i]])
        content    = rawToChar(response$content)
        content_df = read_tsv(content)

        append     = TRUE
        col_names  = FALSE

        if(i == 1){
            append    = FALSE
            col_names = TRUE
        }

        write_tsv(content_df, file = path_out,
                  append = append, col_names = col_names)
    }
    path_out
}

