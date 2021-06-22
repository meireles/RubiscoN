#' Create a lookup table with the accepted species name using `taxize`
#'
#' MOBOTs tropicos sometimes returns a family name. Those get set to NA
#'
#' @param x vector of species names
#' @param chunk max names to submit to TNRS server at once.
#' @param sleep sleep time in seconds before the request for the next chunk is sent
#' @return data frame with the submitted species name and the accepted name
get_accepted_plant_name = function(x, chunk = 1000, sleep = 0){

    families        = unique(c(taxize::apg_families$family,
                               taxize::apg_families$synonym,
                               taxize::apg_families$accepted_name))

    x_capitalized   = str_to_sentence(x)
    x_unique        = unique(x_capitalized)
    x_match_back    = match(x_capitalized, x_unique)

    message("*** ", length(x_unique), " unique names of out ", length(x), " names provided ***")

    parsed          = rgnparser::gn_parse_tidy(x_unique)[["canonicalfull"]]
    parsing_failed  = is.na(parsed)

    parsed[ parsing_failed ] = x_unique[ parsing_failed ]

    parsed = gsub(pattern = "×", replacement = "x", x = parsed)

    split_factor  = sort(rep_len(seq(ceiling(length(parsed)) / chunk), length(parsed)))
    sp_names_list = split(parsed, split_factor)


    matched_list  = lapply(sp_names_list, function(y){

        Sys.sleep(sleep)

        dup       = which(duplicated(y))
        z         = y
        z[dup]    = dup
        otol      = taxize::tol_resolve(z, context_name = "Vascular plants")
        accepted  = otol$unique_name
        unmatched = is.na(accepted)

        if( any(unmatched) ){
            message("unmatched: ", sum(unmatched), "/", chunk)

            mobot_id = 165 # from taxize::gnr_datasources()
            mobot    = gnr_resolve(sci             = z[unmatched],
                                   best_match_only = TRUE,
                                   data_source_ids = mobot_id,
                                   with_context    = TRUE,
                                   canonical       = FALSE)

            m_mobot           = match(mobot$user_supplied_name, z)
            accepted[m_mobot] = mobot$matched_name
        } else {
            message("all names matched")
        }

        if(length(dup) > 0){
            m_dup         = sapply(y[dup], function(l){ which(z == l)  })
            accepted[dup] = accepted[m_dup]
        }

        accepted = rgnparser::gn_parse_tidy(accepted)[["canonicalfull"]]

        accepted[ accepted %in% families ] = NA

        accepted

    })

    accepted_name = do.call(c, matched_list)

    data.frame(original_name  = x,
               parsed_name    = parsed[x_match_back],
               parsing_failed = parsing_failed[x_match_back],
               accepted_name  = accepted_name[x_match_back])
}
