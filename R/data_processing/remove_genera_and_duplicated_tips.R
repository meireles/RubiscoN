
which_genera_and_duplicated_species = function(y){
    which(duplicated(y) | ! base::grepl(y, pattern = " ", fixed = TRUE))
}

remove_genera_and_duplicated_tips = function(x){
    y = x$tip.label
    m = which_genera_and_duplicated_species(y)
    ape::drop.tip(x, m, collapse.singles = TRUE, trim.internal = TRUE)
}

