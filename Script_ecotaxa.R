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
              "vegan", "morphr", "purrr", "imager", "ggrepel", "cowplot", "Nmisc")

# Install and load packages if not already installed
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

## Data Import and Setup ####

# Set the working directory (user should replace with their own path)
setwd("/your/working/directory/path")  # Update with the correct path

# Reads the Ecotaxa file
data <- read_csv("your/ecotaxa/file.csv") 
head(data)

# Reads the Hydrography file
env_data <- read_csv("your/hydrography/file.csv") 
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
env_data <- env_data %>%
  rename(
    Temp = `avg. t`,
    Sal = `avg sal`,
    Fluor = `avg fluor`,
    O2 = `avg o2`,
    Depth_max = `Max Depth`,
    Depth_min = `Min Depth`
  ) %>%
  # Add 'n' prefix to Net column to match the format in data_p
  mutate(Net = paste0("n", Net)) 

# Merge environmental data with processed data
data_p <- data %>%
  left_join(env_data, by = c("cruise" = "Cruise", "moc" = "Tow", "net" = "Net")) %>%
  # Remove 'object_' prefix from column names for simplicity
  rename_with(~ gsub("^object_", "", .x)) 

# Standardize environmental variables
data_p[, c(170:173, 175)] <- decostand(data_p[, c(170:173, 175)], method = "standardize")

# Transform and select morphological variables
spe <- data_p[, 25:91] # Specify the columns of the morphological variables 
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

## Image Processing ####

# Add image paths for each data entry
images_directory <- "imgs/"
data_t1 <- data_t1 %>% mutate(img_path = str_c(images_directory, id, ".jpg"))

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
data_t1$img_path %>% purrr::walk(~ img_chop(.x, bottom = 20))

# Add cropped image paths
cropped_images_directory <- "cropped_imgs/"
data_t1 <- data_t1 %>% mutate(cr_img_path = str_c(cropped_images_directory, id, ".jpg"))

## Principal Components Analysis (PCA) ####

# Prepare data for PCA
weights <- as.data.frame(data_t1[, "density"])
weights <- as.numeric(weights$density)
data_for_pca <- as.data.frame(data_t1[, c("depth_max", "Temp", "Sal", 
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
                       data_t1$cr_img_path, 
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
sil_kmeans <- fviz_nbclust(data_t1[, variables], kmeans, method = "silhouette") 
ggsave("optimal_clusters_silhouette.png", plot = sil_kmeans)
print(sil_kmeans)

# Perform k-means clustering
set.seed(123) # For reproducibility
kmeans_result <- kmeans(data_t1[, variables], centers = 7, nstart = 25, iter.max = 100) 
# Add cluster assignment to data
data_t1$Cluster <- kmeans_result$cluster 

# PCA plot with clusters
var <- res.pca$eig[c(1, 2), 2]
labels <- paste0("PC", c(1, 2), " (", round(var, 1), "%)")
labels <- as.data.frame(labels)
cols <- c("1" = "#543005", "2" = "#8C510A", "3" = "#BF812D", "4" = "#DFC27D", 
          "5" = "#35978F", "6" = "#01665E", "7" = "#003C30")

c1 <- ggplot(data_t1, aes(x = res.pca$ind$coord[, 1], y = res.pca$ind$coord[, 2])) +
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
  ggtitle("PCA with Clusters")

print(c1)
ggsave("pca_clusters_plot.png", plot = c1)

## Density Plots ####

# Day and Night density plots
day_plot <- ggplot(data_t1 %>% filter(moc == "m7"), 
                   aes(x = net, y = density, fill = as.factor(Cluster))) +
  geom_col(position = "stack") +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x = "Depth (m)", y = "Density (ind m³)", fill = "Cluster") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_reverse() +
  scale_x_discrete(
    breaks = unique(data_t1$net), 
    labels = c("1000", "800", "600", "400", "200", "100", "50", "25")) +
  ggtitle("DAY")
ggsave("day_density_plot.png", plot = day_plot)

night_plot <- ggplot(data_t1 %>% filter(moc == "m8"), 
                     aes(x = net, y = density, fill = as.factor(Cluster))) +
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
    breaks = unique(data_t1$net), 
    labels = c("1000", "800", "600", "400", "200", "100", "50", "25")) +
  ggtitle("NIGHT")
ggsave("night_density_plot.png", plot = night_plot)

# Combine day and night plots side by side
combined_plot <- plot_grid(day_plot, night_plot, align = "h", ncol = 2)
print(combined_plot)
ggsave("combined_day_night_density_plot.png", plot = combined_plot)

# Calculate density difference between day and night
density_difference <- data_t1 %>%
  filter(moc %in% c("m7", "m8")) %>%
  group_by(Cluster, net) %>%
  summarise(day_density = sum(density[moc == "m7"]),
            night_density = sum(density[moc == "m8"]),
            density_diff = day_density - night_density) %>%
  ungroup()

# Plot density difference
density_diff_plot <- ggplot(density_difference, 
                            aes(x = net, y = density_diff, fill = as.factor(Cluster))) +
  geom_col(position = "stack") +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x = "Depth (m)", y = "Density Difference (Day - Night) (ind m³)", 
       fill = "Cluster") +
  theme_classic() +
  scale_y_reverse() +
  scale_x_discrete(
    breaks = unique(density_difference$net), 
    labels = c("1000", "800", "600", "400", "200", "100", "50", "25")) +
  ggtitle("Difference in Density (Day - Night)")
print(density_diff_plot)
ggsave("density_difference_plot.png", plot = density_diff_plot)

## Cluster-Specific Density/Abundance Plots ####

# Plot density per cluster in different format 
dens_plot_cl <- ggplot(data_t1, aes(x = net, y = density, 
                                    group = interaction(as.factor(Cluster), moc),  
                                    colour = as.factor(Cluster), linetype = moc, 
                                    fill = as.factor(Cluster), alpha = moc)) +
  stat_summary(geom = "line", fun = sum, size = 1.5) +
  stat_summary(geom = "area", fun = sum) +
  facet_wrap(~Cluster, scales = "free_y", nrow = 1) +
  scale_x_discrete(limits = rev,
                   breaks = unique(data_t1$net), 
                   labels = c( "1000", "800", "600", "400", "200", "100", "50", "25"), 
                   name = "Depth (m)") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_linetype_manual(values = c("solid", "dashed"),
                        breaks = c("m7", "m8"), 
                        labels = c("Day", "Night"), 
                        name = "Day/night") +
  scale_alpha_manual(values = c(0.1, 0.3), 
                     name = "Day/night", 
                     breaks = c("m7", "m8"), 
                     labels = c("Day", "Night")) +
  theme_classic() +
  theme(text = element_text(size = 25, face = "bold"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 3),
        legend.position = "none", 
        axis.text.y = element_text(angle = 90, hjust = 1),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.background = element_rect(fill = "transparent"),
        legend.key.width = unit(3, "line")) +
  ylab(expression(bold("Density ind" ~ m^-3))) +
  guides(alpha = guide_legend(nrow = 1), 
         linetype = guide_legend(nrow = 1))

print(dens_plot_cl)
ggsave("abund_plot_cal.png", width = 40, height = 10, bg = "white")

## Plot Transparency vs. Size with Marginal Distribution ####

# Calculate summary statistics
data_t1$esd_mm <- data_p$esd * 0.0053
data_t1$Mean1 <- data_p$mean
data_t1$Cluster <- as.factor(data_t1$Cluster)

df.mm <- data_t1 %>%
  group_by(Cluster) %>%
  summarise(
    MN.mean = mean(Mean1),
    MN.esd = mean(esd_mm),
    five_mean = quantile(Mean1, probs = 0.05, na.rm = TRUE),
    nintyfive_mean = quantile(Mean1, probs = 0.95, na.rm = TRUE),
    five_major = quantile(esd_mm, probs = 0.05, na.rm = TRUE),
    nintyfive_major = quantile(esd_mm, probs = 0.95, na.rm = TRUE)
  )

# Main plot: Transparency vs. Size with error bars
hs <- ggplot(df.mm, aes(x = MN.esd, y = MN.mean)) +
  geom_point(aes(color = Cluster), size = 5) +
  geom_errorbar(aes(ymin = five_mean, ymax = nintyfive_mean, color = Cluster), width = 0.2, size = 2) +
  geom_errorbarh(aes(xmin = five_major, xmax = nintyfive_major, color = Cluster), height = 6, size = 2) +
  scale_color_manual(values = cols) +
  ylab("Transparency") +
  xlab("Size class (mm)") +
  theme_minimal() +
  theme(text = element_text(size = 20), 
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 21, b = 0, l = 0))) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 5), breaks = c(2, 4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(120, 225), breaks = c(150, 175, 200))

# Marginal distribution for transparency (Mean1)
hs.mean1 <- ggplot(data_t1, aes(x = Mean1, y = density, 
                                group = Cluster, color = Cluster, fill = Cluster)) +
  stat_summary_bin(geom = "area", fun = sum, bins = 10, alpha = 0.3, size = 1) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  ylab(expression(bold("Density ind" ~ m^-3))) +
  clean_theme() + rremove("legend") +
  scale_x_continuous(expand = c(0, 0), limits = c(120, 225), breaks = c(150, 175, 200)) +
  coord_flip()

# Marginal distribution for size (esd_mm)
data_t1$log_size <- log10(data_t1$esd_mm)
hs.esd1 <- ggplot(data_t1, aes(y = density, x = esd_mm, 
                               group = Cluster, color = Cluster, fill = Cluster)) +
  stat_summary_bin(geom = "area", fun = sum, bins = 10, alpha = 0.3, size = 1) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  ylab(expression(bold("Density ind" ~ m^-3))) +
  clean_theme() + rremove("legend") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 6))

# Combine main plot and marginal distributions
p1 <- insert_xaxis_grob(hs, hs.esd1, grid::unit(0.3, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, hs.mean1, grid::unit(0.3, "null"), position = "right")
final_plot_cal <- ggdraw(p2)

print(final_plot_cal)
ggsave("final_plot_cal.png", width = 5, height = 5, bg = "white")
