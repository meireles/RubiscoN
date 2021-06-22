#' Update species names
#'
#' @param x raw data
#' @param sp_field field with species name to be updated
#' @param lookup taxonomic lookup from get_accepted_plant_name()
#'
#' @return data.frame
update_names = function(x, sp_field, lookup) {
    x[[sp_field]] = lookup$accepted_name
    x
}
