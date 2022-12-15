#loading packages
require(ggplot2)
require(geiger)
require(phylotools)
require(viridis)
require(lme4)
require(r2glmm)
require(car)
require(gridExtra)
require(ape)
require(nlme)
if(!require(MCMCglmm)) install.packages("MCMCglmm")
if(!require(caper)) install.packages("caper")
require(MCMCglmm)
require(caper)
require(mulTree)
require(phytools)


############### MCMCglmm ###############
# loading dadi results file
popsize<-read.csv("1_dadi_summary_sp_level.csv")

# z-transformation of "elevation" and "disturbance" in the dataset
z_elevation<-scale(popsize$mean.elev)
z_disturbance<-scale(popsize$mean.disturbance2)
popsize[,"z_elevation"]<-z_elevation
popsize[,"z_disturbance"]<-z_disturbance


### loading fijian ant beast tree
tree_all<-read.nexus("2_Fijian_ants_beast.tre")


mulTree_data <- as.mulTree(data = popsize, tree = tree_all,
                           taxa = "sp_beast2")

### regression formula for all species.
### note. status is a categorical varible. using factor function here.
formula_a=log_nc_na~z_elevation+z_disturbance+factor(status2)+z_disturbance*factor(status2)

mul_priors <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 0.002)))

# Number of interations
nitt <- 240000
# Length of burnin
burnin <- 40000
# Amount of thinning
thin <- 100
mul_parameters <- c(nitt, thin, burnin)

mulTree(mulTree.data = mulTree_data, formula = formula_a, priors = mul_priors,
        parameters = mul_parameters, output = "fijian_ant_mcmc", ESS = 1000,
        chains = 2)
read.mulTree("fijian_ant_mcmc", convergence = TRUE)
all_models <- read.mulTree("fijian_ant_mcmc")
str(all_models)
summarised_results <- summary(all_models, use.hdr = FALSE, cent.tend = mean,
                              prob = c(75, 25))
summarised_results<-summary(all_models)
summarised_results



pdf(file="Fijian_all_mcmc_summary.pdf",height = 5,width = 5)
plot(summarised_results, horizontal = TRUE, ylab = "", cex.coeff = 0.8,
     main = "Posterior distributions", ylim = c(-3,3), cex.terms = 0.5,
     terms = c("Intercept", "elevation", "disturbance", "endemic","exotic","disturbance*endemic","disturbance*exotic","Phylogeny", "Residuals"),
     col = "grey", cex.main = 0.8)
dev.off()


all_mcmc<-data.frame("elevation"=all_models$z_elevation,"disturbance"=all_models$z_disturbance,"endemic"=all_models$`factor(status2)2_endemic`,"exotic"=all_models$`factor(status2)3_Exotic`,"disturbance*endemic"=all_models$`z_disturbance:factor(status2)2_endemic`,"disturbance*exotic"=all_models$`z_disturbance:factor(status2)3_Exotic`,"Phylogeny"=all_models$phylogenetic.variance)

# plot Posterior distributions for each variables
pdf(file="all_mcmc_elevation.pdf",height = 5,width = 5)
ggplot(all_mcmc,aes(x=elevation))+geom_density(fill="grey")+xlim(-3,3)
dev.off()
pdf(file="all_mcmc_disturbance.pdf",height = 5,width = 5)
ggplot(all_mcmc,aes(x=disturbance))+geom_density(fill="grey")+xlim(-3,3)
dev.off()
pdf(file="all_mcmc_endemic.pdf",height = 5,width = 5)
ggplot(all_mcmc,aes(x=endemic))+geom_density(fill="grey")+xlim(-3,3)
dev.off()
pdf(file="all_mcmc_exotic.pdf",height = 5,width = 5)
ggplot(all_mcmc,aes(x=exotic))+geom_density(fill="grey")+xlim(-3,3)
dev.off()
pdf(file="all_mcmc_elevation*exotic.pdf",height = 5,width = 5)
ggplot(all_mcmc,aes(x=elevation.exotic))+geom_density(fill="grey")+xlim(-2,3)
dev.off()
pdf(file="all_mcmc_elevation*wp.pdf",height = 5,width = 5)
ggplot(all_mcmc,aes(x=elevation.wp))+geom_density(fill="grey")+xlim(-2,3)
dev.off()
pdf(file="all_mcmc_disturbance*endemic.pdf",height = 5,width = 5)
ggplot(all_mcmc,aes(x=disturbance.endemic))+geom_density(fill="grey")+xlim(-3,3)
dev.off()
pdf(file="all_mcmc_disturbance*exotic.pdf",height = 5,width = 5)
ggplot(all_mcmc,aes(x=disturbance.exotic))+geom_density(fill="grey")+xlim(-3,3)
dev.off()
pdf(file="all_mcmc_phylogeny.pdf",height = 5,width = 5)
ggplot(all_mcmc,aes(x=Phylogeny))+geom_density(fill="grey")+xlim(-3,3)
dev.off()


################################################

### for endemic species only
popszie_endemic<-popsize[popsize$status=="endemic",]

### mcmcglm with 100 trees
mulTree_data_endemic <- as.mulTree(data = popszie_endemic, tree = tree_all,
                           taxa = "sp_beast2")

formula_endemic=log_nc_na~z_elevation+z_disturbance

mul_priors <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 0.002)))

# Number of interations
nitt <- 240000
# Length of burnin
burnin <- 40000
# Amount of thinning
thin <- 100
mul_parameters <- c(nitt, thin, burnin)

mulTree(mulTree.data = mulTree_data_endemic, formula = formula_endemic, priors = mul_priors,
        parameters = mul_parameters, output = "fijian_ant_endemic_mcmc", ESS = 1000,
        chains = 2)

all_models_endemic <- read.mulTree("fijian_ant_endemic_mcmc")
str(all_models_endemic)
summarised_results_endemic <- summary(all_models_endemic, use.hdr = FALSE, cent.tend = mean,
                              prob = c(75, 25))
summarised_results_endemic<-summary(all_models_endemic)
summarised_results_endemic
pdf(file="endemic_mcmc_summary.pdf",height = 5,width = 5)
plot(summarised_results_endemic, horizontal = TRUE, ylab = "", cex.coeff = 0.8,
     main = "Posterior distributions", ylim = c(-2,2), cex.terms = 0.5,
     terms = c("Intercept", "elevation", "disturbance","Phylogeny", "Residuals"),
     col = "grey", cex.main = 0.8)
dev.off()

endemic_mcmc<-data.frame("elevation"=all_models_endemic$z_elevation,"disturbance"=all_models_endemic$z_disturbance,"Phylogeny"=all_models_endemic$phylogenetic.variance)
pdf(file="endemic_mcmc_elevation.pdf",height = 5,width = 5)
ggplot(endemic_mcmc,aes(x=elevation))+geom_density(fill="grey")+xlim(-2,2)
dev.off()
pdf(file="endemic_mcmc_disturbance.pdf",height = 5,width = 5)
ggplot(endemic_mcmc,aes(x=disturbance))+geom_density(fill="grey")+xlim(-2,2)
dev.off()
pdf(file="endemic_mcmc_Phylogeny.pdf",height = 5,width = 5)
ggplot(endemic_mcmc,aes(x=Phylogeny))+geom_density(fill="grey")+xlim(-2,2)
dev.off()
