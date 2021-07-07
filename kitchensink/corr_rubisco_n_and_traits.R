library("targets")
library("corrplot")
library("phytools")
library("Rphylopars")
library("bayou")

tar_load(rubisco_aa_composition)
tar_load(rubisco_element_composition)
tar_load(try_combined_traits_by_sp)
tar_load(phy_sp_only)

# Intersect species and trim datasets

sp_rub_intersect     = intersect(rubisco_aa_composition$organism,
                                 rubisco_element_composition$organism)
sp_rub_phy_intersect = intersect(sp_rub_intersect, phy_sp_only$tip.label)

sp_intersect         = intersect(sp_rub_phy_intersect,
                                 try_combined_traits_by_sp$AccSpeciesName)

m_aa  = match(sp_intersect, rubisco_aa_composition$organism)
m_ele = match(sp_intersect, rubisco_element_composition$organism)
m_try = match(sp_intersect, try_combined_traits_by_sp$AccSpeciesName)

aa  = rubisco_aa_composition[ m_aa, ]
ele = rubisco_element_composition[ m_ele, ]
try = try_combined_traits_by_sp[ m_try, ]

rownames(aa)  = aa$organism
rownames(ele) = ele$organism
rownames(try) = try$organism

phy = drop.tip(phy_sp_only, setdiff(phy_sp_only$tip.label, sp_intersect))

# Correlation

mat_aa_try = cbind(aa[ , - which(colnames(aa) %in% c("entry_code", "organism" )) ],
                   try[ , c("LfArea_Lflt", "SLA", "N_Dry", "P_Dry")])

mat_ele_try = cbind(ele[ , - which(colnames(ele) %in% c("entry_code", "organism" )) ],
                    try[ , c("LfArea_Lflt", "SLA", "N_Dry", "P_Dry")])

mat_aa_compl  = mat_aa_try[complete.cases(mat_aa_try) , ]
cmat_aa       = cor(mat_aa_compl)
cmat_aa_p     = cor.mtest(mat_aa_compl)

corrplot(corr = cmat_aa, type = "upper", p.mat = cmat_aa_p$p,
         sig.level = 0.005, insig = "blank")

mat_ele_compl = mat_ele_try[complete.cases(mat_ele_try) , ]
cmat_ele      = cor(mat_ele_compl)
cmat_ele_p    = cor.mtest(mat_ele_compl)

corrplot(corr = cmat_ele, type = "upper", p.mat = cmat_ele_p$p,
         sig.level = 0.005, insig = "blank")

# Phylo Stuff
x   = setNames( mat_ele_compl$N, rownames(mat_ele_compl) )
#x   = log(x)
#x   = x[ x > mean(x) - 4 * sd(x) & x < mean(x) + 4 * sd(x) ]

p = drop.tip(phy, setdiff(phy$tip.label, names(x)))
p = multi2di(p, random = FALSE)
p = force.ultrametric(p)

p$edge.length[p$edge.length < 1e-9] = 1e-9

fit = Rphylopars::anc.recon(x, p)

pdf("kitchensink/aaa.pdf")
contMap(p, x, method = "user", type = "fan", lwd = c(0.5, 2), anc.states = fit,
        fsize = c(0.05, 1), outline = FALSE)
dev.off()


########################################
# Bayou
########################################


prior = make.prior(p,
                   dists = list(dalpha = "dhalfcauchy",
                                dsig2  = "dhalfcauchy",
                                dk     = "cdpois",
                                dtheta = "dnorm"),
                   param = list(dalpha = list(scale = 0.1),
                                dsig2  = list(scale = 0.1),
                                dk     = list(lambda = 10, kmax = 50),
                                dsb    = list(bmax = 1, prob = 1),
                                dtheta = list(mean = mean(x), sd = 1.5 * sd(x)))
                   )

prior_sim = priorSim(prior, p, plot = TRUE)
startpars = prior_sim$pars[[1]]

mcmc_ou = bayou.makeMCMC(tree = p, dat = x, prior = prior, model = "OU",
                         startpar = startpars, outname = "model_ou", plot.freq = 10000)

mcmc_ou$run(500000)

chain_ou = mcmc_ou$load()
summary(chain_ou)
plot(chain_ou, auto.layout = TRUE)


par(mfrow = c(1, 3))
plotSimmap.mcmc(chain_ou, burnin = 0.7, pp.cutoff = 0.3, cex = 0.1)
plotBranchHeatMap(p, chain_ou, "theta", burnin = 0.7, pal = cm.colors, cex = 0.1)
phenogram.density(p, x, burnin = 0.7, chain_ou, pp.cutoff = 0.3, cex = 0.1)




