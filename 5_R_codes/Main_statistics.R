library(tidyverse)
library(brms)
library(MCMCglmm)
library(tidybayes)
mycolors <- c("chartreuse3", "dodgerblue1", "gold")
#
#
#
# import species-level ecological data
grouped <- read_csv("./data/fsc_species_level_summary.csv") %>% 
  select(sp_beast2, status, disturbance_sp, elevation_sp) %>%
  mutate(disturbance_sp = 5 - disturbance_sp) 

# import data on population size changes, m1=honeybee mutation rate
dat <- read_csv("data/fsc_population.csv") %>% mutate(log_nc_na = log10(NCUR_m1/NANC_m1)) %>%
  left_join(grouped)

dat$status <- factor(dat$status, levels=c("endemic", "widespread Pacific native", 'exotic'))  

# average population size changes to species level
dat_grouped <- dat %>% group_by(sp_beast2) %>% summarise(log_nc_na = mean(log_nc_na), t_ave = mean(T_m1), missing=mean(mean_missing), depth=mean(Average_Depth),disturbance_sp = disturbance_sp[1], elevation_sp = elevation_sp[1], status = status[1] )

dat_grouped %>% ggplot(aes(disturbance_sp, log_nc_na, color = status)) + geom_point() + stat_smooth(method = "lm")
#
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

phylo <- calTree(ape::read.tree(tree_files[1]), dat$sp_beast2)
A <- ape::vcv.phylo(phylo)

#Full model, including the sequencing depth and missing data to access the effect of allelic dropout
model <- brm(
  log_nc_na ~ (depth+missing+disturbance_sp + elevation_sp) * status  + (1|gr(sp_beast2, cov = A)) ,
  # t_ave ~ (depth+missing+disturbance_sp + elevation_sp) * status  + (1|gr(sp_beast2, cov = A)), # model for t
  data = dat_grouped,
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
  phylo <- calTree(ape::read.tree(file), dat$sp_beast2)
  A <- ape::vcv.phylo(phylo)
  model <- update(model, newdata = dat_grouped, data2 = list(A = A), recompile = FALSE)
}
# save model
saveRDS(model, "model_m1.rds") 
#
#
#
model <- readRDS("model_m1.rds")
#model_t<-readRDS("model_m1_t.rds")
#
#
#
summary(model, digits = 3)

# update model to remove the sequence depth and missing for the main results.
model2 <- update(model, . ~ .  -depth - depth:status -missing - missing:status)
summary(model2)

hyp <- "sd_sp_beast2__Intercept^2 / (sd_sp_beast2__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(model2, hyp, class = NULL))
bayes_R2(model2)
#
#
# Plot figure 3
ce <- conditional_effects(model2, "disturbance_sp:status")

ce2 <- conditional_effects(model2, "status", conditions = list(disturbance_sp = 0)) # compute marginal effects for species in the primary forest

p_main <- 
  plot(ce, plot = FALSE, points = T)[[1]] + scale_fill_manual(values = mycolors) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = status), alpha = 0.01, colour = NA) +
  labs(x = "Species average disturbance score", y = expression(paste("Population size change (", log[10], frac("present","past"),")"))) +
  scale_color_manual(values=mycolors) +
  theme_bw() +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  guides(color = "none") +
  theme(legend.position = c(.8, .2), legend.title = element_blank(), legend.background = element_rect(colour="grey", size=0.5, linetype="solid") ) +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  ylim(-5.5, 5.5) 
p_main
# X-axis density plot (on top)
p_x_density <- ggplot(dat_grouped, aes(x = disturbance_sp, fill = status)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = mycolors) +
  xlim(0,3)+
  theme_void() +  # Remove axis elements for cleaner visualization
  theme(legend.position = "none")

# Y-axis density plot (on the right, **forced to be narrower**)
p_y_density <- ggplot(dat_grouped, aes(y = log_nc_na, fill = status)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = mycolors) +
  ylim(-5.5,5.5)+
  theme_void() +
  theme(legend.position = "none")

# Empty spacer for correct alignment
empty_plot <- ggplot() + theme_void()

# Combine plots with correct proportions
final_plot <- (p_x_density + plot_spacer()) / 
  (p_main + p_y_density) +
  plot_layout(heights = c(0.5, 4), widths = c(4, 0.5))
final_plot
ggsave("figure3.pdf", height=5, width=7)


# Plot posterior distributions
posterior <- as_draws_df(model)  # Extract posterior samples
posterior_long <- posterior %>%
  pivot_longer(cols = starts_with("b_") | starts_with("sd_"), names_to = "Parameter", values_to = "Estimate")
#rename messy parameter names to more readable labels.
posterior_long <- posterior_long %>%
  mutate(Parameter = case_when(
    Parameter == "b_Intercept" ~ "Intercept (Endemic, Disturbance = 0)",
    Parameter == "b_disturbance_sp" ~ "Effect of Disturbance",
    Parameter == "b_elevation_sp" ~ "Effect of Elevation",
    Parameter == "b_statuswidespreadPacificnative" ~ "Widespread Pacific Native (vs. Endemic)",
    Parameter == "b_statusexotic" ~ "Exotic Species (vs. Endemic)",
    Parameter == "b_disturbance_sp:statuswidespreadPacificnative" ~ "Disturbance × Widespread Pacific Native",
    Parameter == "b_disturbance_sp:statusexotic" ~ "Disturbance × Exotic Species",
    Parameter == "b_elevation_sp:statuswidespreadPacificnative" ~ "Elevation × Widespread Pacific Native",
    Parameter == "b_elevation_sp:statusexotic" ~ "Elevation × Exotic Species",
    Parameter == "sd_sp_beast2__Intercept" ~ "Phylogenetic Variance (Species-level)"
  )) 
#Reorder Parameters (Top to Bottom)
posterior_long <- posterior_long %>%
  mutate(Parameter = factor(Parameter, levels = c(
    "Phylogenetic Variance (Species-level)",  # Random effect variance
    "Elevation × Exotic Species", 
    "Elevation × Widespread Pacific Native", 
    "Effect of Elevation", 
    "Disturbance × Exotic Species", 
    "Disturbance × Widespread Pacific Native", 
    "Exotic Species (vs. Endemic)", 
    "Widespread Pacific Native (vs. Endemic)", 
    "Effect of Disturbance", 
    "Intercept (Endemic, Disturbance = 0)"
  )))
# plot
ggplot(posterior_long, aes(x = Estimate, y = Parameter, fill = Parameter)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) +
  scale_fill_viridis_d() +
  labs(title = "Posterior Distributions of Model Parameters",
       x = "Posterior distributions") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("posterior_distributions2.pdf", height=7, width=7)



