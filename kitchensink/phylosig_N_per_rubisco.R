library("targets")
library("phytools")
library("Rphylopars")
library("bayou")

tar_load(final_data)

clade = "All"
data  = final_data[[clade]]

normalize_by_seq_length = TRUE


phy   = data$tree



try   = data$try_data[ , c("LfArea_Lflt", "SLA", "N_Dry", "P_Dry") ] # "VcM_Dry"
ele   = data$elem_comp[ , c("C", "N")]
aa    = data$aa_comp
aa    = aa[ , - which(colnames(aa) %in% c("entry_code", "organism") )]

# Normalize by sequence length?
if(normalize_by_seq_length){
    s   = rowSums(aa)
    aa  = aa / s
    ele = ele / s
}

ele   = cbind(ele, "C:N" = ele$C / ele$N )

x     = setNames(ele$C, data$elem_comp$organism)

fit   = Rphylopars::anc.recon(x, phy)

pdf("kitchensink/phylo_elem.pdf")
contMap(phy, x, method = "user", type = "fan", lwd = c(0.5, 2), anc.states = fit,
        fsize = c(0.05, 1), outline = FALSE)
dev.off()


########################################
# Bayou
########################################


prior = make.prior(phy,
                   dists = list(dalpha = "dhalfcauchy",
                                dsig2  = "dhalfcauchy",
                                dk     = "cdpois",
                                dtheta = "dnorm"),
                   param = list(dalpha = list(scale = 0.005), #list(scale = 0.1),
                                dsig2  = list(scale = 0.005), #list(scale = 0.1),
                                dk     = list(lambda = 5, kmax = 20), #list(lambda = 10, kmax = 50),
                                dsb    = list(bmax = 1, prob = 1),
                                dtheta = list(mean = mean(x), sd = 1.5 * sd(x)))
)

prior_sim = priorSim(prior, phy, plot = TRUE)
startpars = prior_sim$pars[[1]]

mcmc_ou = bayou.makeMCMC(tree = phy, dat = x, prior = prior, model = "OU",
                         startpar = startpars, outname = "model_ou", plot.freq = 10000)

mcmc_ou$run(100000)

chain_ou = mcmc_ou$load()
summary(chain_ou)
plot(chain_ou, auto.layout = TRUE)


par(mfrow = c(1, 3))
plotSimmap.mcmc(chain_ou, burnin = 0.5, pp.cutoff = 0.5, cex = 0.1)
plotBranchHeatMap(phy, chain_ou, "theta", burnin = 0.5, pal = cm.colors, cex = 0.1)
phenogram.density(phy, x, burnin = 0.5, chain_ou, pp.cutoff = 0.5, cex = 0.1)




