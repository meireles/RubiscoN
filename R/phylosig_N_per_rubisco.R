# tar_load(rubisco_element_composition)
# tar_load(phy_sp_only)
#
# X = aggregate(N ~ organism, data = rubisco_element_composition, mean)
# x = setNames(X$N, X$organism)
#
# i = intersect(phy_sp_only$tip.label, names(x))
# x = x[i]
#
# phy = drop.tip(phy_sp_only, setdiff(phy_sp_only$tip.label, i))
#
