# ## Trait data
#
# Here is a list requested traits and their TRY IDs
#
# 1. Leaf Area: 1, 3108, 3109, 3110, 3111, 3112, 3113, 3114
# 2. Specific Leaf Area: 11, 3115, 3116, 3117, 3085, 3086
# 3. Leaf lifespan: 12
# 4. Leaf nitrogen (N) content per leaf area: 50
# 5. Leaf nitrogen (N) content per leaf dry mass: 14
# 6. Leaf photosynthesis pathway: 22
# 7. Leaf photosynthesis carboxylation capacity (Vcmax) per leaf area (Farquhar model): 186
# 8. Leaf photosynthesis carboxylation capacity (Vcmax) per leaf dry mass (Farquhar model): 185
# 9. Leaf phosphorus (P) content per leaf area: 51
# 10. Leaf phosphorus (P) content per leaf dry mass: 15

aggregate_traits_by_sp = function(x){

    nice_mean  = try_keep_txt(mean)
    nice_merge = function(a, b){
        merge(a, b, all = TRUE)
    }

    x = x[ , c("AccSpeciesName", "ObservationID", "TraitID", "TraitName",
               "ValueKindName", "StdValue", "UnitName")]

    # Remove genera and NA species
    x = x[ ! is.na(x$AccSpeciesName) & grepl(" ", x$AccSpeciesName, fixed = TRUE), ]

    # 22 is all NAs
    y = aggregate(StdValue ~ AccSpeciesName + TraitName + TraitID,
                  data = x, FUN = nice_mean, na.rm = TRUE)

    # 22 is all NAs and thus removed
    # 3086 is not included in 'l' and thus removed
    l = c("LfArea_Uncert" = 3114, "LfArea_WholeLf" = 3108, "LfArea_WholeLf" = 3110,
          "LfArea_WholeLf" = 3112, "LfArea_Lflt" = 3109, "LfArea_Lflt" = 3111,
          "LfArea_Lflt" = 3113, "SLA" = 3115, "SLA" = 3116, "SLA" = 3117,
          "Lf_Life" = 12, "N_Area" = 50,  "N_Dry" = 14, "P_Area" = 51, "P_Dry" = 15,
          "VcM_Area" = 186, "VcM_Dry" = 185)

    i = setNames(unique(names(l)), unique(names(l)))

    z = lapply(i, function(n){
        message(n)
        a = y[ y$TraitID %in% l[ names(l) == n ] ,  ]
        b = aggregate(StdValue ~ AccSpeciesName, data = a, FUN = nice_mean, na.rm = TRUE)
        names(b)[2] = n
        b
    })

    m = Reduce(nice_merge, z)
    m
}



# combine_try_data_by_specimen = function(x){
#
#     x = x[ , c("AccSpeciesName", "ObservationID", "TraitID", "TraitName",
#                "ValueKindName", "StdValue", "UnitName")]
#
#     x = x[ ! is.na(x$AccSpeciesName), ]
#
#     p = aggregate(AccSpeciesName ~ ObservationID, x, unique)
#     r = unique(x$ObservationID)
#     s = unique(x$TraitName)
#
#     y = data.frame(matrix(nrow     = length(r),
#                           ncol     = length(s),
#                           dimnames = list(NULL, s)),
#                    check.names = FALSE)
#
#     y = cbind(p, y)
#
#     for(i in s){
#         message(i)
#         z       = x[ x$TraitName == i ,  ]
#         m       = match(z$ObservationID, y$ObservationID)
#         y[m, i] = z$StdValue
#     }
#
#     y
# }
