#### R script for analysis of paired D/N MOCNESS tows, after scanning in zooscan and uploaded to EcoTaxa. 
#### Data uses the ecotaxa csv export file and the hydrographic summary file from the MOCNESS software (BESS or SIO; they should be adapted to the example file). 
#### Script developed by Rocio Rodriguez Perez, rbrodri3@asu.edu

## General Notes: 
# This script is designed to analyze differences in community composition between day and night within a specific group. 
# The current setup focuses on data from a single tow at a specific location, with "m7" representing daytime data and "m8" representing nighttime data. 
# It is important to adjust these parameters based on your specific needs, including the tow (moc) identifiers, the number of clusters, color schemes, and data handling procedures.
# The selection of morphological variables should be revisited and adjusted when applying the script to a different dataset. The correlation matrix may yield different results, highlighting different highly associated variables. Careful consideration should be given to the selection of these variables to ensure accurate analysis. 
# This script is flexible and can be modified accordingly to suit different datasets or analysis requirements.

## Load and install required packages ####

# List of required packages
packages <- c("tidyr", "dplyr", "readr", "stringr", "data.table", "ggplot2", 
              "readxl", "ggpubr", "gridExtra", "RColorBrewer", "colorspace", 
              "FactoMineR", "factoextra", "gginnards", "cellWise", "corrplot", 
              "vegan", "morphr", "purrr", "imager", "ggrepel", "cowplot", "Nmisc", "grid")

# Install and load packages if not already installed
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

## Data Import and Setup ####

# Set the working directory (user should replace with their own path)
setwd("/Users/rociorodriguez/ASU Dropbox/Rocio Rodriguez Perez/ZoopGroup_LAJ/Projects/Gradients_20-8162/CIREN/Rocio")  # Update with the correct path

# Reads the Ecotaxa file
data <- read_csv("ecotaxa_export_5421_20241205_1741.csv") 
head(data)

# Reads the Hydrography file
env_data <- read_csv("Gradients_MOCNESS_net_hydrography.csv") 
head(env_data)

## Data Handling ####

# Note: Modify the column names and keys in steps below according to the format of your data
# Process data:  
data <- data %>%
  # Split 'object_id' into multiple columns
  separate(object_id, c("cruise", "moc", "net", "fraction"), extra = "drop", remove = FALSE) %>% 
  # Calculate density and abundance.m2
  mutate(density = acq_sub_part / sample_tot_vol,
         abundance.m2 = (acq_sub_part / sample_tot_vol) * 
           (object_depth_max - object_depth_min))

# Rename columns in environmental data (for simplicity) 
env_data1 <- env_data %>%
  rename(
    Temp = `temp`,
    Sal = `sal`,
    Fluor = `avg fluor`,
    O2 = `avg o2`,
    Depth_max = `max depth`,
    Depth_min = `min depth`,
    D.N = `D/N`
  ) %>%
  # Add 'n' prefix to Net column to match the format in data_p
  mutate(Cruise = tolower(Cruise),
         net = paste0("n", net), 
         tow = paste0("m", tow)) 



# Merge environmental data with processed data
data_p <- data %>%
  left_join(env_data1, by = c("cruise" = "Cruise", "moc" = "tow", "net" = "net")) %>%
  # Remove 'object_' prefix from column names for simplicity
  rename_with(~ gsub("^object_", "", .x)) 

# Standardize environmental variables
data_p[, c("Temp", "Sal", "Fluor", "O2", "Depth_max", "Depth_min")] <- decostand(data_p[, c("Temp", "Sal", "Fluor", "O2", "Depth_max", "Depth_min")], method = "standardize")

# Transform and select morphological variables
spe <- data_p[, 23:89] # Specify the columns of the morphological variables 
YJ_trans <- transfo(spe, type = "YJ")$Y # Perform a Yeo-Johnson transformation on the morphological variables 

data_t1 <- data_p %>%
  # Remove the original morphological variables from the data
  select(-one_of(colnames(spe))) %>%
  # Bind the transformed morphological variables to the remaining data
  bind_cols(YJ_trans)

# Store the transformed data in a separate object
YJ_data <- YJ_trans

## Morphological Variable Selection ####

# Computes the correlation matrix
cor_matrix <- cor(YJ_data) 

# Plots the correlation matrix
corrplot(cor_matrix, method = "color", type = "upper", number.cex = 0.4) 

# Define morphological variables for clustering and PCA
# Highly correlated variables (|r| > 0.90) were removed to avoid redundancy
variables <- c("mean", "mode", "min", "max", "angle", "circ.", 
               "kurt", "fractal", "nb1", "symetrievc", "fcons", 
               "thickr", "esd", "elongation", "meanpos", "sr", 
               "feretareaexc", "perimferet", "cdexc")

data_t1 <- data_t1 %>%
  mutate(Region = case_when(
    Station %in% c(7, 11) ~ "Subarctic",
    Station %in% c(15, 21, 23) ~ "Central",
    Station %in% c(27, 34) ~ "Subtropical",
    TRUE ~ NA_character_
  ))
# Subset the data to only include rows where the "annotation_category" column is equal to "Calanoida"
Cal_data <- subset(data_t1, annotation_category == "Calanoida")

## Image Processing ####

# Add image paths for each data entry
images_directory <- "imgs/"
Cal_data <- Cal_data %>% mutate(img_path = str_c(images_directory, id, ".jpg"))

# Create directory for cropped images
if (!dir.exists("cropped_imgs")) {
  dir.create("cropped_imgs")
}

# Function to crop images
img_chop <- function(img_path, bottom = 0, save_path = "cropped_imgs") {
  # Load the image
  img <- load.image(img_path) 
  # Get the width of the image
  w <- width(img) 
  # Get the height of the image
  h <- height(img) 
  # Crop the image at the bottom
  cropped_img <- img[1:w, 1:(h-bottom), , , drop = FALSE] 
  # Get the image name
  img_name <- basename(img_path)
  # Construct the save path
  save_img_path <- file.path(save_path, img_name) 
  # Save the cropped image
  save.image(cropped_img, save_img_path) 
}

# Crop images
Cal_data$img_path %>% purrr::walk(~ img_chop(.x, bottom = 20))

# Add cropped image paths
cropped_images_directory <- "cropped_imgs/"
Cal_data <- Cal_data %>% mutate(cr_img_path = str_c(cropped_images_directory, id, ".jpg"))

## Principal Components Analysis (PCA) ####

# Prepare data for PCA
weights <- as.data.frame(Cal_data[, "density"])
weights <- as.numeric(weights$density)
data_for_pca <- as.data.frame(Cal_data[, c("depth_max", "Temp", "Sal", 
                                          "O2", "Fluor", variables)])

# Convert all columns in data_for_pca to numeric
data_for_pca <- data_for_pca %>%
  mutate(across(everything(), as.numeric))

# Perform PCA
res.pca <- PCA(data_for_pca, 
               scale.unit = TRUE, 
               quanti.sup = 1:5, 
               graph = FALSE, 
               row.w = weights)

# Get PCA variables and create PCA plot with images
pca.vars <- rbind(res.pca$var$coord, res.pca$quanti.sup$coord) %>% as.data.frame
images <- ggmorph_tile(res.pca, 
                       Cal_data$cr_img_path, 
                       dimensions = c(1, 2), 
                       steps = 9, n_imgs = 10, 
                       scale = 0.003, 
                       adjust_grey = TRUE)

# PCA plot with arrows and labels
images1 <- images +
  geom_segment(data = pca.vars, aes(x = 0, y = 0, xend = Dim.1 * 4.5, yend = Dim.2 * 4.5), 
               # Add arrows to the plot
               arrow = arrow(length = unit(1/2, 'picas'))) + 
  # Add text labels
  geom_text_repel(data = pca.vars, aes(x = Dim.1 * 5, y = Dim.2 * 5, label = rownames(pca.vars))) + 
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) 

print(images1)
# Save PCA plot
ggsave("pca_plot_with_images.png", plot = images1)

## Clustering ####

# Determine optimal number of clusters using the silhouette method
sil_kmeans <- fviz_nbclust(Cal_data[, variables], kmeans, method = "silhouette") 
ggsave("optimal_clusters_silhouette.png", plot = sil_kmeans)
print(sil_kmeans)

# Perform k-means clustering
set.seed(123) # For reproducibility
kmeans_result <- kmeans(Cal_data[, variables], centers = 6, nstart = 25, iter.max = 100) 
# Add cluster assignment to data
Cal_data$Cluster <- kmeans_result$cluster 

# PCA plot with clusters
var <- res.pca$eig[c(1, 2), 2]
labels <- paste0("PC", c(1, 2), " (", round(var, 1), "%)")
labels <- as.data.frame(labels)
cols <- c("1" = "#543005", "2" = "#8C510A", "3" = "#BF812D", "4" = "#DFC27D", 
          "5" = "#35978F",  "6" = "#003C30")

c1 <- ggplot(Cal_data, aes(x = res.pca$ind$coord[, 1], y = res.pca$ind$coord[, 2])) +
  geom_point(aes(color = as.factor(Cluster))) + 
  theme_classic() +
  xlab(as.character(labels[1, ])) + ylab(as.character(labels[2, ])) +
  theme(text = element_text(size = 30, face = "bold"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        legend.position = "bottom", 
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 10), nrow = 1, byrow = TRUE)) +
  scale_color_manual(values = cols, name = "Clusters") +
  ggtitle("PCA with Clusters")+
  facet_wrap(~Region)

print(c1)
ggsave("pca_clusters_plot.png", plot = c1)

#### 10. Generate Abundance Plots for Each Region ####
# Loop through each region
for (region in unique(Cal_data$Region)) {
  
  # Filter the data by region
  region_data <- Cal_data %>% filter(Region == region)
  
  # Summarize the data to get the total density per 'moc' and then calculate the average
  summarized_data <- region_data %>%
    group_by(moc) %>%
    summarize(total_density = sum(density, na.rm = TRUE)) %>%
    ungroup() %>%
    summarize(avg_density = mean(total_density, na.rm = TRUE))
  
  # Print or use the average density as needed
  print(paste("Average density for region", region, ":", summarized_data$avg_density))
  
  # Create the day plot for the region
  day_plot <- ggplot(region_data %>% filter(D.N == "D"), aes(x = net, y = density, fill = as.factor(Cluster))) +
    geom_col(position = "stack") +
    scale_fill_manual(values = cols) +
    coord_flip() +
    labs(x = "Depth (m)", y = "Density (ind m³)", fill = "Cluster") +
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_reverse() +
    scale_x_discrete(
      breaks = unique(region_data$net), 
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25")) +
    ggtitle(paste0("DAY - ", region))
  
  # Create the night plot for the region
  night_plot <- ggplot(region_data %>% filter(D.N == "N"), aes(x = net, y = density, fill = as.factor(Cluster))) +
    geom_col(position = "stack") +
    scale_fill_manual(values = cols) +
    coord_flip() +
    labs(x = "Depth (m)", y = "Density (ind m³)", fill = "Cluster") +
    theme_classic() +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.line.y = element_blank(), 
          plot.title = element_text(hjust = 1)) +
    scale_x_discrete(
      breaks = unique(region_data$net), 
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25")) +
    ggtitle(paste0("NIGHT - ", region))
  
  # Combine the day and night plots side by side
  combined_plot <- plot_grid(day_plot, night_plot, align = "h", ncol = 2)
  
  # Save the combined plot for the region
  ggsave(paste0("combined_day_night_density_plot_", region, ".png"), plot = combined_plot)
}


# Calculate density difference between day and night
density_difference <- Cal_data %>%
  filter(D.N %in% c("D", "N")) %>%
  group_by(Cluster, net) %>%
  summarise(day_density = sum(density[D.N == "D"]),
            night_density = sum(density[D.N == "N"]),
            density_diff = day_density - night_density) %>%
  ungroup()

# Loop through each region to create Density Difference plots
for (region in unique(Cal_data$Region)) {
  
  # Filter the data by region
  region_data <- Cal_data %>% filter(Region == region)
  
  # Calculate density difference (Day - Night)
  density_difference <- region_data %>%
    group_by(net, Cluster) %>%
    summarize(
      day_density = sum(density[D.N == "D"], na.rm = TRUE),
      night_density = sum(density[D.N == "N"], na.rm = TRUE)
    ) %>%
    mutate(density_diff = day_density - night_density)
  
  # Plot density difference
  density_diff_plot <- ggplot(density_difference, aes(x = net, y = density_diff, fill = as.factor(Cluster))) +
    geom_col(position = "stack") +
    scale_fill_manual(values = cols) +
    coord_flip() +
    labs(x = "Depth (m)", y = "Density Difference (Day - Night) (ind m³)", fill = "Cluster") +
    theme_classic() +
    scale_y_reverse() +
    scale_x_discrete(
      breaks = unique(density_difference$net), 
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25")) +
    ggtitle(paste0("Difference in Density (Day - Night) - ", region))
  
  # Save density difference plot for the region
  ggsave(paste0("density_difference_plot_", region, ".png"), plot = density_diff_plot, width = 10, height = 6)
}

# Loop through each region to create Density per Cluster plots
for (region in unique(Cal_data$Region)) {
  
  # Filter the data by region
  region_data <- Cal_data %>% filter(Region == region)
  
  # Plot density per cluster
  dens_plot_cl <- ggplot(region_data, aes(x = net, y = density, group = interaction(as.factor(Cluster), D.N), 
                                          colour = as.factor(Cluster), linetype = D.N, fill = as.factor(Cluster), alpha = D.N)) +
    stat_summary(geom = "line", fun = sum, size = 1.5) +
    stat_summary(geom = "area", fun = sum) +
    facet_wrap(~Cluster, scales = "free_y", nrow = 1) +
    scale_x_discrete(
      limits = rev,
      breaks = unique(region_data$net), 
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25"), 
      name = "Depth (m)") +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_linetype_manual(values = c("solid", "dashed"),
                          breaks = c("D", "N"), 
                          labels = c("Day", "Night"), 
                          name = "Day/night") +
    scale_alpha_manual(values = c(0.1, 0.3), 
                       name = "Day/night", 
                       breaks = c("D", "N"), 
                       labels = c("Day", "Night")) +
    theme_classic() +
    theme(text = element_text(size = 25, face = "bold"), 
          panel.border = element_rect(colour = "black", fill = NA, size = 3),
          legend.position = "none", 
          axis.text.y = element_text(angle = 90, hjust = 1),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x = element_text(angle = 180, vjust = 1),
          plot.background = element_rect(fill = "transparent"),
          legend.key.width = unit(3, "line")) +
    ylab(expression(bold("Density ind" ~ m^-3))) +
    guides(alpha = guide_legend(nrow = 1), 
           linetype = guide_legend(nrow = 1))
  
  # Save density per cluster plot for the region
  ggsave(paste0("abund_plot_cal_", region, ".png"), plot = dens_plot_cl, width = 40, height = 10, bg = "white")
}

####Plot transparency ~ size with marginal distribution#####
Cal_data$esd_mm<- data_p$esd*.0053
Cal_data$Mean1<- data_p$mean
Cal_data$Cluster<-as.factor(Cal_data$Cluster)

# Initialize a list to store the plots
plots_list <- list()

# Loop through each region
for (region in unique(Cal_data$Region)) {
  
  # Filter the data by region
  region_data <- Cal_data %>% filter(Region == region)
  
  # Summarize the data for the region
  df.mm <- region_data %>%
    group_by(Cluster) %>%
    summarise(
      MN.mean = mean(Mean1),
      MN.esd = mean(esd_mm),
      five_mean = quantile(Mean1, probs = 0.05, na.rm = TRUE),
      nintyfive_mean = quantile(Mean1, probs = 0.95, na.rm = TRUE),
      five_major = quantile(esd_mm, probs = 0.05, na.rm = TRUE),
      nintyfive_major = quantile(esd_mm, probs = 0.95, na.rm = TRUE)
    )
  
  # Create hs plot for the region
  hs <- ggplot(data=df.mm, aes(x=MN.esd, y=MN.mean)) +
    geom_point(aes(color=Cluster), size=5) +
    geom_errorbar(aes(ymin=five_mean, ymax=nintyfive_mean, color=Cluster), width=0.2, size=2) +
    geom_errorbarh(aes(xmin=five_major, xmax=nintyfive_major, color=Cluster), height=6, size=2) +
    scale_color_manual(values=cols) +
    ylab("Transparency") +
    xlab("Size class (mm)") +
    theme_minimal() +
    theme(text = element_text(size=20), 
          legend.position = "none",
          axis.title.y = element_text(margin=margin(t=0, r=21, b=0, l=0)))
  
  # Create hs.mean1 plot for the region
  hs.mean1 <- ggplot(region_data, aes(x=Mean1, y=density, group=Cluster, color=Cluster, fill=Cluster)) +
    stat_summary_bin(geom="area", fun=sum, bins = 10, alpha=0.3, size=1) +
    scale_color_manual(values=cols) +
    scale_fill_manual(values=cols) +
    theme_classic() +
    ylab(expression(bold("Density ind" ~ m^-3))) +
    clean_theme() + 
    rremove("legend") +
    coord_flip()
  
  # Create hs.esd1 plot for the region
  hs.esd1 <- ggplot(region_data, aes(y=density, x=esd_mm, group=Cluster, color=Cluster, fill=Cluster)) +
    stat_summary_bin(geom="area", fun=sum, bins = 10, alpha=0.3, size=1) +
    scale_color_manual(values=cols) +
    scale_fill_manual(values=cols) +
    theme_classic() +
    ylab(expression(bold("Density ind" ~ m^-3))) +
    clean_theme() + 
    rremove("legend")
  
  # Insert hs.esd1 as the top x-axis and hs.mean1 as the right y-axis
  p1 <- insert_xaxis_grob(hs, hs.esd1, grid::unit(.3, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, hs.mean1, grid::unit(.3, "null"), position = "right")
  
  # Final combined plot for the region
  final_plot_cal <- ggdraw(p2)
  
  # Save the plot in the list
  plots_list[[region]] <- final_plot_cal
  
  # Save the final plot for the region as an image
  ggsave(paste0("final_plot_cal_", region, ".png"), plot = final_plot_cal, width = 10, height = 10, bg = "white")
}

# Display the list of plots
plots_list
