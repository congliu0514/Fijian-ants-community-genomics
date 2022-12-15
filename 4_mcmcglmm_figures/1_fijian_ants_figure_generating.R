#loading packages
require(ggplot2)
require(viridis)
### Loading datasets
species_data<-read.csv("0_Fijian_ants_species_level_data_checked2208.csv",stringsAsFactors = FALSE)
pop_fst<-read.csv("1_Fijian_ants_population_fst_data.csv")
dadi_all<-read.csv("2_dadi_summary_sp_level.csv",stringsAsFactors = F)

########################################################################################################################

# 1. Ploting figure 2a: Species ecology
# 1.1 scatter plot
pdf(file="fig_2a_species_ecology.pdf",width = 5, height = 5)
ggplot(data=species_data,aes(x=-1*species_data$mean.disturbance,y=species_data$mean.elev,colour=as.factor(species_data$endemic),shape=as.factor(species_data$endemic)))+geom_point(size=3)+scale_color_manual(values=c("#336699", "#339999", "#CC6666"))+theme(legend.position = "none")
dev.off()
# 1.2 distribution of elevation
pdf(file="fig_2a_elevation_distribution.pdf",width = 5, height = 5)
ggplot(species_data,aes(x=species_data$mean.elev,fill=as.factor(species_data$endemic)))+geom_density(alpha=0.4)
dev.off()
# 1.3 distribution of disturbance
pdf(file="fig_2a_disturbance_distribution.pdf",width = 5, height = 5)
ggplot(species_data,aes(x=-1*species_data$mean.disturbance,fill=as.factor(species_data$endemic)))+geom_density(alpha=0.4)
dev.off()
######################################################################################################################## 

# 2. Ploting figure 2b: FST~ distance
# 2.1 scatter plot
pdf(file="fig_2b_fst_geo-distance.pdf",width = 5, height = 5)
ggplot(data=pop_fst,aes(x=pop_fst$distance,y=pop_fst$fst, colour=as.factor(pop_fst$endemic),shape=as.factor(pop_fst$endemic)))+geom_point(size=3)+scale_color_manual(values=c("#336699", "#339999", "#CC6666"))+theme(legend.position = "none")+geom_smooth(method=lm,se=T)
dev.off()
# 2.2 distribution of FST.
pdf(file="fig_2b_fst_distribution.pdf",width = 5, height = 5)
ggplot(pop_fst,aes(x=pop_fst$fst,fill=as.factor(pop_fst$endemic)))+geom_density(alpha=0.4)
dev.off()
########################################################################################################################

# 3. Ploting figure 2C: log_nc_na ~ disturbance
# 3.1 scatter plot
pdf(file="fig_2c_population_size_disturbance_3.pdf",width = 5, height = 5)
ggplot(data=dadi_all,aes(x=-1*dadi_all$mean.disturbance,y=dadi_all$log_nc_na,colour=as.factor(dadi_all$status2),shape=as.factor(dadi_all$status2)))+geom_point(size=3)+scale_color_manual(values=c("#336699", "#339999", "#CC6666"))+theme(legend.position = "none")+geom_smooth(method=lm,se=T)
dev.off()
# 3.2 distribution of log_nc_na
pdf(file="fig_2c_population_size_distribution.pdf",width = 5, height = 5)
ggplot(dadi_all,aes(x=log_nc_na,fill=as.factor(dadi_all$status2)))+geom_density(alpha=0.4)
dev.off()
########################################################################################################################

# 4. Ploting figure 1b: distributions of t.
# 4.1 for endemic species
dadi_End<-dadi_all[dadi_all$status=="endemic",]
pdf(file="density_time_endemic.pdf",width = 5,height = 5)
number_of_cases=nrow(dadi_End) # plot the distributions of population contraction/expansion together relative to the total number of endemic species
ggplot(dadi_End,aes(x=t_years,y=(..count..)/number_of_cases*1000,fill=as.factor(population.decline)))+geom_density(alpha=0.4)+xlim(0,3000)
dev.off()

# 4.2 for widespread pacific native species
dadi_WP<-dadi_all[dadi_all$status=="widespread pacific native",]
pdf(file="density_time_WP2.pdf",width = 5,height = 5)
ggplot(dadi_WP,aes(x=t_years))+geom_density(fill="grey")+xlim(0,3000)
dev.off() 

# 4.3 for exotic species
dadi_Exotic<-dadi_all[dadi_all$status=="Exotic",]
pdf(file="density_time_exotic.pdf",width = 5,height = 5)
ggplot(dadi_Exotic,aes(x=t_years))+geom_density(fill="grey")+xlim(0,3000)
dev.off()
