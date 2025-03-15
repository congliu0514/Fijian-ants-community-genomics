library(tidyverse)
library(brms)
library(MCMCglmm)
library(tidybayes)
mycolors <- c("chartreuse3", "dodgerblue1", "gold")
#
#
#
# import species-level ecological data
data1 <- read_csv("./data/fsc_species_level_summary.csv") %>% 
  select(sp_beast2, status, disturbance_sp, elevation_sp) %>%
  mutate(disturbance_sp = 5 - disturbance_sp) 


#
calTree <- function(tree, taxa) { # remove taxa not found in sp_beast2 and calibrate the tree
  tree <- drop.tip(tree, setdiff(tree$tip.label, taxa))
  root_node <- c("Proceratium.relictum_EGP0091C09_Fiji", "Hypoponera.eutrepta_EGP0093B09_Fiji", "Hypoponera.monticola_EGP0094C04_Fiji", "Leptogenys.FJE02_EGP0095E10_Fiji", "Leptogenys.letilae_EGP0095D08_Fiji", "Leptogenys.FJE05_EGP0096A07_Fiji", "Leptogenys.fjn02_EGP0095E11_Fiji", "Leptogenys.fjn01_EGP0095F08_Fiji", "Odontomachus.angulatus_EGP0096C03_Fiji", "Odontomachus.simillimus_EGP0098D03_Fiji", "Anochetus.graeffei_EGP0091H05_Fiji", "Pseudoponera.stigma_EGP0099C08_Fiji")
  tree <- root(tree, outgroup = root_node, resolve.root = TRUE)
  node <- c(
    getMRCA(tree, tip = c("Pheidole.roosevelti_EGP0016D06_Fiji","Rogeria.stigmatica_EGP0035C04_Fiji") ), # Myrmicinae: 51-72 MY
    getMRCA(tree, tip = c("Paraparatrechina.oceanica_EGP0138A07_Fiji","Colobopsis.vitiensis_EGP0060A03_Fiji") ), #Formicinae: 51-71 MY
    getMRCA(tree, tip = c("Hypoponera.eutrepta_EGP0093B09_Fiji","Pseudoponera.stigma_EGP0099C08_Fiji") ) #  Ponerinae: 61-84 MY
  )
  age.min <- c(51,51,61)
  age.max <- c(72,71,84)
  soft.bounds <- c(FALSE,FALSE,FALSE)
  mycalibration <- data.frame(node, age.min, age.max, soft.bounds)
  tree <- chronos(tree, lambda = 1, model = "relaxed", calibration = mycalibration, control = chronos.control() )
  return(tree)
}
#
#
#
tree_files <- list.files("data/100bootstrap_Fijian_all_specimens", pattern = "Fijian_ants_bootstrap_Tree_\\d{2}\\.tre", full.names = TRUE)

phylo <- calTree(ape::read.tree(tree_files[1]), data1$sp_beast2)
A <- ape::vcv.phylo(phylo)
model <- brm(
  disturbance_sp ~ status  + (1|gr(sp_beast2, cov = A)) ,
  data = data1,
  data2 = list(A = A),
  family = gaussian(),
  prior = c(
    prior(normal(0,10), "b"),
    prior(normal(0,50), "Intercept"),
    prior(student_t(3,0,20), "sd"),
    prior(student_t(3,0,20), "sigma")
  ),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  sample_prior = TRUE, chains = 2, cores = 2, 
  iter = 4000, warmup = 2000
)

for (file in tree_files[-1]) {
  phylo <- calTree(ape::read.tree(file), data1$sp_beast2)
  A <- ape::vcv.phylo(phylo)
  model <- update(model, newdata = data1, data2 = list(A = A), recompile = FALSE)
}
#
#
#
saveRDS(model, "model_figure1_disturbance.rds")
summary(model)



# Extract posterior samples from the model
posterior_samples <- as_draws_df(model)

# Compute the difference between exotic and widespread Pacific native
posterior_samples <- posterior_samples %>%
  mutate(diff_exotic_native = b_statusexotic - b_statuswidespreadPacificnative)

# Summarize the posterior difference
posterior_summary <- posterior_samples %>%
  summarize(
    mean = mean(diff_exotic_native),
    lower_95 = quantile(diff_exotic_native, 0.025),
    upper_95 = quantile(diff_exotic_native, 0.975)
  )

# Print the result
print(posterior_summary)
#
ggplot(posterior_samples, aes(x = diff_exotic_native)) +
  geom_density(fill = "dodgerblue", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Posterior Difference: Exotic vs. Widespread Pacific Native",
       x = "Difference (statusexotic - statuswidespreadPacificnative)",
       y = "Density") +
  theme_minimal()

