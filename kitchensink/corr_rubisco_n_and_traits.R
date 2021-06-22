library("targets")
library("corrplot")

tar_load(rubisco_element_composition)
tar_load(try_combined_traits_by_sp)

# Intersect species and trim datasets
sp_intersect = intersect(rubisco_element_composition$organism,
                         try_combined_traits_by_sp$AccSpeciesName)

m_rub = match(sp_intersect, rubisco_element_composition$organism)
m_try = match(sp_intersect, try_combined_traits_by_sp$AccSpeciesName)

rub = rubisco_element_composition[ m_rub , ]
try = try_combined_traits_by_sp[ m_try , ]


# Correlation
mat = cbind(rub[ , c("C", "H", "N", "O", "S")],
            try[ , c("SLA", "N_Area", "P_Area", "VcM_Area")])

rownames(mat) = rub$organism

cmat = cor(mat[complete.cases(mat), ])

corrplot(cmat)
