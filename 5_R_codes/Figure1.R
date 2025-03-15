# Load required libraries
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggridges)
library(patchwork)
library(dunn.test)

# Step 1: Read the sorted tree
tree <- read.tree("./Data/fijian_beast_tree2.tre")
tree <- ladderize(tree, right = FALSE)

# Step 2: Read the original dataset
data <- read.csv("./Data/fijian_tree_species_ecology.csv", stringsAsFactors = FALSE)

# Step 3: Read the new dataset (fst_d values)
fst_data <- read.csv("./Data/population_fst_data.csv", stringsAsFactors = FALSE)

# Step 4: Ensure `tree_id` format matches tree tip labels in both datasets
tree$tip.label <- gsub("'", "", tree$tip.label)
data$tree_id <- gsub("_", ".", data$tree_id)  
fst_data$tree_id <- gsub("_", ".", fst_data$tree_id)  
tree$tip.label <- gsub("_", ".", tree$tip.label)  

# Step 5: Merge the new dataset with the original dataset based on `tree_id`
data <- merge(data, fst_data, by = "tree_id", all.x = TRUE)  # Ensure all original species are retained

# Step 6: Convert `endemic_class` to a factor
data$endemic_class <- as.factor(data$endemic_class)

# Step 7: Reshape dataset into long format, adding `fst_d`
tree_data <- data %>%
  select(tree_id, endemic_class, Disturbance, Elevation, fst_d) %>%
  pivot_longer(cols = c(Disturbance, Elevation, fst_d), names_to = "Variable", values_to = "Value")

# Ensure `tree_data` follows `tree$tip.label` order
tree_data$tree_id <- factor(tree_data$tree_id, levels = rev(tree$tip.label))

# Check that tree_data now has data
print(table(tree_data$Variable))  # Debugging step

# Define colors for endemic classes
endemic_colors <- c("widespread pacific native" = "blue",
                    "endemic" = "green",
                    "Exotic" = "yellow")  

# Step 8: Extract the root age (time depth of the tree)
max_time <- max(node.depth.edgelength(tree))

# Step 9: Plot the **time-scaled phylogenetic tree**
p_tree <- ggtree(tree, size = 0.7) + 
  geom_tiplab(align = TRUE, size = 3, hjust = -0.1) +
  theme_tree2() +
  ggtitle("Phylogenetic Tree with Time Scale") +
  scale_x_continuous(name = "Time (millions of years ago)", limits = c(0, max_time))  

# Step 10: **Density plot for Disturbance**
p_disturbance <- ggplot(tree_data %>% filter(Variable == "Disturbance"),
                        aes(x = Value, y = tree_id, fill = endemic_class)) +
  geom_density_ridges(alpha = 0.7, scale = 1) +  
  scale_fill_manual(values = endemic_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 4)) +  
  theme_minimal() +
  labs(x = "Disturbance", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ggtitle("Disturbance Density")

# Step 11: **Density plot for Elevation**
p_elevation <- ggplot(tree_data %>% filter(Variable == "Elevation"),
                      aes(x = Value, y = tree_id, fill = endemic_class)) +
  geom_density_ridges(alpha = 0.7, scale = 1) +  
  scale_fill_manual(values = endemic_colors, guide = "none") +  
  scale_x_continuous(limits = c(0, 1500)) +  
  theme_minimal() +
  labs(x = "Elevation", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +  
  ggtitle("Elevation Density")

# Step 12: **Density plot for Fst_d (Genetic Differentiation)**
p_fst_d <- ggplot(tree_data %>% filter(Variable == "fst_d"),
                  aes(x = Value, y = tree_id, fill = endemic_class)) +
  geom_density_ridges(alpha = 0.7, scale = 1) +  
  scale_fill_manual(values = endemic_colors, guide = "none") +  
  theme_minimal() +xlim(0,0.005)+
  labs(x = "Fst_d", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +  
  ggtitle("Fst_d Density")

# Step 13: Arrange tree and plots in four columns
final_plot <- p_tree | p_disturbance | p_elevation | p_fst_d  

print(final_plot)

# Step 14: Save as PDF
pdf(file="figure1_with_fst_d.pdf", width=25, height=15)  # Increased width for 4 columns
print(final_plot)
dev.off()


#### boxplots for figure 1
library(ggplot2)
library(dplyr)

# Step 1: Read the dataset
data <- read.csv("fijian_tree_species_ecology.csv", stringsAsFactors = FALSE)

data_fst<-read.csv("population_fst_data.csv")
# Step 2: Convert `endemic_class` to a factor **and reverse the order**
data$endemic_class <- factor(data$endemic_class, levels = rev(c("endemic", "widespread pacific native", "Exotic")))
data_fst$endemic_class<- factor(data_fst$endemic_class, levels = rev(c("endemic", "widespread pacific native", "Exotic")))
# Define colors for endemic classes
endemic_colors <- c("endemic" = "green",
                    "widespread pacific native" = "blue",
                    "Exotic" = "yellow")

# Step 3: **Boxplot for Disturbance** with corrected order
p_disturbance_boxplot <- ggplot(data, aes(x = Disturbance, y = endemic_class, fill = endemic_class)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = endemic_colors, guide = "none") +
  theme_minimal()+ 
  labs(x = "Disturbance", y = "Endemic Class") +
  ggtitle("Disturbance Boxplot")

# Step 4: **Boxplot for Elevation** with corrected order
p_elevation_boxplot <- ggplot(data, aes(x = Elevation, y = endemic_class, fill = endemic_class)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = endemic_colors, guide = "none") +
  theme_minimal() +
  labs(x = "Elevation", y = "Endemic Class") +
  ggtitle("Elevation Boxplot")

# Step 4-2: **Boxplot for fst** with corrected order
p_fst_boxplot <- ggplot(data_fst, aes(x = fst_d, y = endemic_class, fill = endemic_class)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = endemic_colors, guide = "none") +
  theme_minimal() +xlim(0,0.005)+
  labs(x = "FST", y = "Endemic Class") +
  ggtitle("FST Boxplot")


# Step 5: Print the corrected boxplots
print(p_disturbance_boxplot)
print(p_elevation_boxplot)
print(p_fst_boxplot)

# Save as PDF
pdf(file="corrected_boxplots_fst.pdf", width=5, height=5)
#print(p_disturbance_boxplot)
#print(p_elevation_boxplot)
print(p_fst_boxplot)
dev.off()







