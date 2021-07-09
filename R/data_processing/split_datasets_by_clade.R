
split_data_by_clade_and_repackage = function(phy, try, aa, ele){

    ########################################
    # Subset tree and datasets
    ########################################

    # To list the labeled nodes on the tree
    # phy$node.label[ ! ( grepl("mrcaott", phy$node.label) | phy$node.label == "") ]


    # Remove misplaced fern everywhere
    phy = ape::drop.tip(phy = phy, tip = "Asplenium rutifolium")

    # Pick node labels to split and find their node number
    clades_split = c("Caryophyllales", "Poales", "Asterales", "Ericales", "Laurales",
                     "Malpighiales", "Fabales", "Sapindales", "Apiales", "Fagales",
                     "Myrtales")

    clades_node  = sapply(clades_split, grep, x = phy$node.label)
    clades_node  = clades_node + ape::Ntip(phy)

    # Add "All" (entire tree) and Gymnosperms (their mrca isn't named) to the mix
    clades_node  = c("All"         = ape::Ntip(phy) + 1L,
                     "Gymnosperms" = ape::getMRCA(phy, c("Pinus taeda", "Ginkgo biloba")),
                     clades_node)

    ########################################
    # Subset tree and datasets
    ########################################

    ####################
    # subset data helper
    ####################

    dat_subset   = function(tree, dat, colname_match){
        dat[ match(tree$tip.label, dat[ , colname_match]),   ]
    }


    # Make a list of subtrees
    subtree_list = sapply(clades_node, ape::extract.clade, phy = phy,
                          simplify = FALSE)


    # Subset datasets
    try_list = sapply(subtree_list, dat_subset, dat = try,
                      colname_match = "AccSpeciesName", simplify = FALSE)

    aa_list = sapply(subtree_list, dat_subset, dat = aa,
                     colname_match = "organism", simplify = FALSE)

    ele_list = sapply(subtree_list, dat_subset, dat = ele,
                      colname_match = "organism", simplify = FALSE)

    ########################################
    # Repackage
    ########################################

    dat = lapply(setNames(names(clades_node), names(clades_node)), function(x){
        list(tree      = subtree_list[[x]],
             try_data  = try_list[[x]],
             aa_comp   = aa_list[[x]],
             elem_comp = ele_list[[x]])
    })

    dat
}
