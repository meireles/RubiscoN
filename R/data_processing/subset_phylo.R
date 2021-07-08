

subset_phylo = function(tree, sp_keep, brlens_tol = 1e-9){
    phy = ape::drop.tip(tree, setdiff(tree$tip.label, sp_keep))
    phy = ape::multi2di(phy, random = FALSE)
    phy = phytools::force.ultrametric(phy)
    phy$edge.length[phy$edge.length < brlens_tol-9] = brlens_tol
    ladderize(phy)
}
