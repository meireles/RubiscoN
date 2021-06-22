remove_genera_and_duplicated_tips = function(x){
    y = x$tip.label
    m = which(duplicated(y) | ! base::grepl(y, pattern = " ", fixed = TRUE))
    ape::drop.tip(x, m, collapse.singles = TRUE, trim.internal = TRUE)
}
