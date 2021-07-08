library("targets")
library("ape")

tar_load(phy_sp_only)
tar_load(try_combined_traits_by_sp)

node_keep = grep("Caryophyllales", phy_sp_only$node.label) + ape::Ntip(phy_sp_only)
tips_keep = ape::extract.clade(phy_sp_only, node_keep)$tip.label
dat       = try_combined_traits_by_sp[try_combined_traits_by_sp$AccSpeciesName %in% tips_keep , ]

write.csv(dat, "kitchensink/try_caryophyllales_all.csv", row.names = FALSE)

