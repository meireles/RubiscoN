# library("phytools")
# library("Rphylopars")
# library("bayou")
#
# tar_load(rubisco_element_composition_intersect)
# tar_load(phy_sp_only_intersect)
#
#
# X = aggregate(N ~ organism, data = rubisco_element_composition, mean)
# x = setNames(X$N, X$organism)
#
# i = intersect(phy_sp_only$tip.label, names(x))
# x = x[i]
#
# phy = drop.tip(phy_sp_only, setdiff(phy_sp_only$tip.label, i))
#
#
#
#
# x   = setNames( mat_ele_compl$N, rownames(mat_ele_compl) )
# #x   = log(x)
# #x   = x[ x > mean(x) - 4 * sd(x) & x < mean(x) + 4 * sd(x) ]
#
# p = drop.tip(phy, setdiff(phy$tip.label, names(x)))
# p = multi2di(p, random = FALSE)
# p = force.ultrametric(p)
#
# p$edge.length[p$edge.length < 1e-9] = 1e-9
#
# fit = Rphylopars::anc.recon(x, p)
#
# pdf("kitchensink/aaa.pdf")
# contMap(p, x, method = "user", type = "fan", lwd = c(0.5, 2), anc.states = fit,
#         fsize = c(0.05, 1), outline = FALSE)
# dev.off()
#
#
# ########################################
# # Bayou
# ########################################
#
#
# prior = make.prior(p,
#                    dists = list(dalpha = "dhalfcauchy",
#                                 dsig2  = "dhalfcauchy",
#                                 dk     = "cdpois",
#                                 dtheta = "dnorm"),
#                    param = list(dalpha = list(scale = 0.1),
#                                 dsig2  = list(scale = 0.1),
#                                 dk     = list(lambda = 10, kmax = 50),
#                                 dsb    = list(bmax = 1, prob = 1),
#                                 dtheta = list(mean = mean(x), sd = 1.5 * sd(x)))
# )
#
# prior_sim = priorSim(prior, p, plot = TRUE)
# startpars = prior_sim$pars[[1]]
#
# mcmc_ou = bayou.makeMCMC(tree = p, dat = x, prior = prior, model = "OU",
#                          startpar = startpars, outname = "model_ou", plot.freq = 10000)
#
# mcmc_ou$run(500000)
#
# chain_ou = mcmc_ou$load()
# summary(chain_ou)
# plot(chain_ou, auto.layout = TRUE)
#
#
# par(mfrow = c(1, 3))
# plotSimmap.mcmc(chain_ou, burnin = 0.7, pp.cutoff = 0.3, cex = 0.1)
# plotBranchHeatMap(p, chain_ou, "theta", burnin = 0.7, pal = cm.colors, cex = 0.1)
# phenogram.density(p, x, burnin = 0.7, chain_ou, pp.cutoff = 0.3, cex = 0.1)
#
#
#
#
