#' Read and trim TRY data
#'
#' @param path_in path for TRY's zip file
#'
#' @return tibble
read_and_trim_raw_try_data = function(path_in){

    o = dirname(path_in)

    # Unzip data
    zip::unzip(path_in, exdir = o)

    # get a handle for the data file
    p = gsub(".zip", ".txt", path_in)

    # Read data
    x = vroom::vroom(p)

    # Delete the data file
    file.remove(p)

    # Remove metadata
    # That is, everything that does not have a trait id.
    x = x[ ! is.na(x$TraitID), ]

    # Remove useless columms
    keep = setdiff(colnames(x), c("LastName", "FirstName", "DatasetID",
                                  "Dataset", "SpeciesName", "DataID",
                                  "OriglName", "Comment", "...28"))
    x = x[ , keep]

    # remove special chars from the species names
    # see https://stackoverflow.com/questions/9934856/removing-non-ascii-characters-from-data-files
    x$AccSpeciesName = iconv(x$AccSpeciesName, "latin1", "ASCII", sub = "")

    x
}
