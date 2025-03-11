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

# note for unix systems: ensure that the `cmake` system library is installed
# note that `morphr` is not on CRAN: https://github.com/jiho/morphr

packages <- c(
  "tidyr", "dplyr", "readr", "stringr", "data.table", "ggplot2",
  "readxl", "ggpubr", "gridExtra", "RColorBrewer", "colorspace",
  "FactoMineR", "factoextra", "gginnards", "cellWise", "corrplot",
  "vegan", "morphr", "purrr", "imager", "ggrepel", "cowplot", "Nmisc", "grid"
)

# Install and load packages if not already installed

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

## Data Import and Setup ####

# Set the working directory (user should replace with their own path)
# setwd("/Users/rociorodriguez/ASU Dropbox/Rocio Rodriguez Perez/ZoopGroup_LAJ/Projects/Gradients_20-8162/CIREN/Rocio")  # Update with the correct path

env <- "sol"
env <- "local"

if (env == "sol") {
  this_path = "/scratch/srearl/ciren/rocio"
} else if (env == "local") {
  this_path = here::here("rocio")
}

# Reads the Ecotaxa file
# data <- readr::read_csv(here::here("rocio", "ecotaxa_export_5421_20241205_1741.csv"))
# data <- readr::read_csv("/scratch/srearl/ciren/rocio/ecotaxa_export_5421_20241205_1741.csv")
data <- readr::read_csv(here::here(this_path, "ecotaxa_export_5421_20241205_1741.csv"))

# head(data)

# Reads the Hydrography file
# env_data <- readr::read_csv("/scratch/srearl/ciren/rocio/Gradients_MOCNESS_net_hydrography.csv")
env_data <- readr::read_csv(here::here(this_path, "Gradients_MOCNESS_net_hydrography.csv"))
# head(env_data)


## Data Handling ####

# Note: Modify the column names and keys in steps below according to the format of your data
# Process data:
data <- data |>
  # Split 'object_id' into multiple columns
  tidyr::separate(
    col    = object_id,
    into   = c("cruise", "moc", "net", "fraction"),
    extra  = "drop",
    remove = FALSE
  ) |>
  # Calculate density and abundance.m2
  dplyr::mutate(
    density      = acq_sub_part / sample_tot_vol,
    abundance.m2 = (acq_sub_part / sample_tot_vol) * (object_depth_max - object_depth_min)
  )

# Rename columns in environmental data (for simplicity)
env_data1 <- env_data |>
  dplyr::rename(
    Temp      = `temp`,
    Sal       = `sal`,
    Fluor     = `avg fluor`,
    O2        = `avg o2`,
    Depth_max = `max depth`,
    Depth_min = `min depth`,
    D.N       = `D/N`
  ) |>
  # Add 'n' prefix to Net column to match the format in data_p
  dplyr::mutate(
    Cruise = tolower(Cruise),
    net    = paste0("n", net),
    tow    = paste0("m", tow)
  )


# Merge environmental data with processed data
data_p <- data |>
  dplyr::left_join(
    y  = env_data1,
    by = c(
      "cruise" = "Cruise",
      "moc"    = "tow",
      "net"    = "net"
    )
  ) |>
  # Remove 'object_' prefix from column names for simplicity
  dplyr::rename_with(~ gsub("^object_", "", .x))

# Standardize environmental variables
data_p[, c("Temp", "Sal", "Fluor", "O2", "Depth_max", "Depth_min")] <- vegan::decostand(
  x      = data_p[, c("Temp", "Sal", "Fluor", "O2", "Depth_max", "Depth_min")],
  method = "standardize"
)

# Transform and select morphological variables
spe <- data_p[, 23:89] # Specify the columns of the morphological variables 
YJ_trans <- cellWise::transfo(spe, type = "YJ")$Y # Perform a Yeo-Johnson transformation on the morphological variables 

# SRE: cellWise produces a list of which the Y object is the data matrix

data_t1 <- data_p |>
  # Remove the original morphological variables from the data
  dplyr::select(-one_of(colnames(spe))) |>
  # Bind the transformed morphological variables to the remaining data
  dplyr::bind_cols(YJ_trans)

# Store the transformed data in a separate object
YJ_data <- YJ_trans

# SRE: at this point we have (1) environmental variables that have been
# standardized with deostand, (2) morphological variables that have been
# transormed per Yeo-Johnson both as (a) part of data_t1 and (b) as a
# standalone object (YJ_data), and (3) associated tow details.

## Morphological Variable Selection ####

# Computes the correlation matrix
cor_matrix <- cor(YJ_data) 

# Plots the correlation matrix
corrplot::corrplot(
  corr       = cor_matrix,
  method     = "color",
  type       = "upper",
  number.cex = 0.4
)

# Define morphological variables for clustering and PCA
# Highly correlated variables (|r| > 0.90) were removed to avoid redundancy
variables <- c(
  "mean",
  "mode",
  "min",
  "max",
  "angle",
  "circ.",
  "kurt",
  "fractal",
  "nb1",
  "symetrievc",
  "fcons",
  "thickr",
  "esd",
  "elongation",
  "meanpos",
  "sr",
  "feretareaexc",
  "perimferet",
  "cdexc"
)

data_t1 <- data_t1 |>
  dplyr::mutate(
    Region = dplyr::case_when(
      Station %in% c(7, 11) ~ "Subarctic",
      Station %in% c(15, 21, 23) ~ "Central",
      Station %in% c(27, 34) ~ "Subtropical",
      TRUE ~ NA_character_
    )
  )

# Subset the data to only include rows where the "annotation_category" column is equal to "Calanoida"
Cal_data <- subset(data_t1, annotation_category == "Calanoida")

# Image Processing -------------------------------------------------------------

# Add image paths for each data entry

if (env == "sol") {
  images_directory <- here::here(this_path, "imgs")
} else if (env == "local") {
  images_directory <- "~/Desktop/images/imgs/"
}

# images_directory <- "/scratch/srearl/imgs/"
# images_directory <- "imgs/"

Cal_data <- Cal_data |>
  dplyr::mutate(img_path = stringr::str_c(images_directory, id, ".jpg"))

# SRE:save to parquet
# arrow::write_parquet(
#   x = Cal_data,
#   sink = "cal_data.parquet"
# )

# --- start sbatch image_chop

from_parquet <- arrow::read_parquet("cal_data.parquet")
identical(Cal_data, from_parquet)

# Create directory for cropped images
if (!dir.exists("cropped_imgs")) {
  dir.create("cropped_imgs")
}

# Function to crop images
img_chop <- function(img_path, bottom = 0, save_path = "cropped_imgs") {
  
  # Load the image
  img <- imager::load.image(img_path) 
  
  # Get the width of the image
  w <- imager::width(img) 
  
  # Get the height of the image
  h <- imager::height(img) 
  
  # Crop the image at the bottom
  cropped_img <- img[1:w, 1:(h-bottom), , , drop = FALSE] 
  # print(cropped_img)
  
  # Get the image name
  img_name <- basename(img_path)
  
  # Construct the save path
  save_img_path <- file.path(save_path, img_name)
  
  # Save the cropped image
  imager::save.image(cropped_img, save_img_path)
  
}

# just one
# img_chop(Cal_data$img_path[1], bottom = 20)

# Crop images
# SRE: this would be a good call to parallelize
Cal_data$img_path |>
  purrr::walk(~ img_chop(.x, bottom = 20))

# --- end sbatch image_chop

# Add cropped image paths
cropped_images_directory <- "~/cropped_imgs/"

Cal_data <- Cal_data |>
  dplyr::mutate(cr_img_path = stringr::str_c(cropped_images_directory, id, ".jpg"))

# SRE:save to parquet
arrow::write_parquet(
  x = Cal_data,
  sink = "/scratch/srearl/cal_data.parquet"
)

from_parquet <- arrow::read_parquet("/scratch/srearl/cal_data.parquet")
identical(Cal_data, from_parquet)

# Principal Components Analysis (PCA) ------------------------------------------

# Prepare data for PCA
weights <- as.data.frame(Cal_data[, "density"])
colnames(weights) <- "density" # SRE: added to ensure as.numeric works on the generated column
weights <- as.numeric(weights$density)


data_for_pca <- as.data.frame(
  Cal_data[, c(
    "depth_max", "Temp", "Sal",
    "O2", "Fluor", variables
  )]
)

# Convert all columns in data_for_pca to numeric
data_for_pca <- data_for_pca |>
  dplyr::mutate(across(everything(), as.numeric)) # |> dplyr::slice(1:10000)

# SRE: PCA on environmental and morphological variables only on samples
# identified as Calanoida. Curious mixing of dependent and independent
# variables in a single analysis. Maybe a question for Rocio.

# SRE:save to parquet
# arrow::write_parquet(
#   x = data_for_pca,
#   sink = "/scratch/srearl/data_for_pca.parquet"
# )

# not sure if as.data.frame is needed other than ensuring that identical is TRUE
# from_parquet_pca <- as.data.frame(
#   x = arrow::read_parquet("/scratch/srearl/data_for_pca.parquet")
# )
# identical(data_for_pca, from_parquet_pca)

# Perform PCA
res.pca <- FactoMineR::PCA(
  X          = data_for_pca,
  scale.unit = TRUE,
  quanti.sup = 1:5,
  graph      = FALSE,
  row.w      = weights
)

# Get PCA variables and create PCA plot with images
pca.vars <- rbind(res.pca$var$coord, res.pca$quanti.sup$coord) |>
  as.data.frame()

# SRE: I can run the PCA but am unable to plot using the supplied code without
# the image data.

# need sbatch to complete!
images <- morphr::ggmorph_tile(
  space       = res.pca,
  imgs        = Cal_data$cr_img_path,
  dimensions  = c(1, 2),
  steps       = 9,
  n_imgs      = 10,
  scale       = 0.003,
  adjust_grey = TRUE
)

# Get PCA variables and create PCA plot with images
pca.vars <- rbind(
  res.pca$var$coord,
  res.pca$quanti.sup$coord
) |>
  as.data.frame()

# PCA plot with arrows and labels
images1 <- images +
  geom_segment(
    data = pca.vars, aes(x = 0, y = 0, xend = Dim.1 * 4.5, yend = Dim.2 * 4.5),
    # Add arrows to the plot
    arrow = arrow(length = unit(1 / 2, "picas"))
  ) +
  # Add text labels
  geom_text_repel(data = pca.vars, aes(x = Dim.1 * 5, y = Dim.2 * 5, label = rownames(pca.vars))) +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

print(images1)
# Save PCA plot
ggsave("pca_plot_with_images.png", plot = images1)


# Clustering -------------------------------------------------------------------

# SRE: pinch point?

# Determine optimal number of clusters using the silhouette method
Rprofmem("Rprofmem.out")
sil_kmeans <- factoextra::fviz_nbclust(
  x          = Cal_data[c(1:1000), variables],
  # x          = Cal_data[, variables],
  FUNcluster = kmeans,
  method     = "silhouette"
)
Rprofmem(NULL)

p <- profmem::profmem({
  
  sil_kmeans <- factoextra::fviz_nbclust(
    x          = Cal_data[c(1:1000), variables],
    # x          = Cal_data[, variables],
    FUNcluster = kmeans,
    method     = "silhouette"
  )
  
})


prof <- profvis::profvis({
  
  sil_kmeans <- factoextra::fviz_nbclust(
    x          = Cal_data[c(1:1000), variables],
    # x          = Cal_data[, variables],
    FUNcluster = kmeans,
    method     = "silhouette"
  )
  
})


# print(p, expr = FALSE)

arrow::write_parquet(
  x = Cal_data[, variables],
  sink = "cal_data_vars.parquet"
)

# another approach but also fails
clus_gap <- cluster::clusGap(
  x = Cal_data[c(1:1000), variables],
  # x = Cal_data[, variables],
  FUN = kmeans,
  K.max = 25
)

gsave("optimal_clusters_silhouette.png", plot = sil_kmeans)

print(sil_kmeans)

# Perform k-means clustering

set.seed(123) # For reproducibility

kmeans_result <- kmeans(
  x        = Cal_data[, variables],
  centers  = 6,
  nstart   = 25,
  iter.max = 100
)

# Add cluster assignment to data
Cal_data$Cluster <- kmeans_result$cluster

# PCA plot with clusters
var    <- res.pca$eig[c(1, 2), 2]
labels <- paste0("PC", c(1, 2), " (", round(var, 1), "%)")
labels <- as.data.frame(labels)
cols   <- c(
  "1" = "#543005",
  "2" = "#8C510A",
  "3" = "#BF812D",
  "4" = "#DFC27D",
  "5" = "#35978F",
  "6" = "#003C30"
)

c1 <- ggplot2::ggplot(
  data    = Cal_data,
  mapping = ggplot2::aes(
    x = res.pca$ind$coord[, 1],
    y = res.pca$ind$coord[, 2]
  )
) +
  ggplot2::geom_point(ggplot2::aes(color = as.factor(Cluster))) +
  ggplot2::theme_classic() +
  ggplot2::xlab(as.character(labels[1, ])) +
  ggplot2::ylab(as.character(labels[2, ])) +
  ggplot2::theme(
    text = ggplot2::element_text(
      size = 30,
      face = "bold"
    ),
    panel.border = ggplot2::element_rect(
      colour = "black",
      fill   = NA,
      size   = 2
    ),
    legend.position = "bottom",
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      override.aes = list(size = 10),
      nrow         = 1,
      byrow        = TRUE
    )
  ) +
  ggplot2::scale_color_manual(
    values = cols,
    name   = "Clusters"
  ) +
  ggplot2::ggtitle("PCA with Clusters") +
  ggplot2::facet_wrap(~Region)

print(c1)
ggsave("pca_clusters_plot.png", plot = c1)


#### 10. Generate Abundance Plots for Each Region ####

# Loop through each region
for (region in unique(Cal_data$Region)) {
  
  # Filter the data by region
  region_data <- Cal_data |>
    dplyr::filter(Region == region)
  
  # Summarize the data to get the total density per 'moc' and then calculate the average
  summarized_data <- region_data |>
    dplyr::group_by(moc) |>
    dplyr::summarize(total_density = sum(density, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::summarize(avg_density = mean(total_density, na.rm = TRUE))
  
  # Print or use the average density as needed
  print(paste("Average density for region", region, ":", summarized_data$avg_density))

  # Create the day plot for the region
  day_plot <- ggplot2::ggplot(
    data = region_data |>
      dplyr::filter(D.N == "D"),
    mapping = ggplot2::aes(
      x = net,
      y = density,
      fill = as.factor(Cluster)
    )
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x    = "Depth (m)",
      y    = "Density (ind m³)",
      fill = "Cluster"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_x_discrete(
      breaks = unique(region_data$net),
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25")
    ) +
    ggplot2::ggtitle(paste0("DAY - ", region))

  # Create the night plot for the region
  night_plot <- ggplot2::ggplot(
    data = region_data |>
      dplyr::filter(D.N == "N"),
    mapping = ggplot2::aes(
      x = net,
      y = density,
      fill = as.factor(Cluster)
    )
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x    = "Depth (m)",
      y    = "Density (ind m³)",
      fill = "Cluster"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.line.y  = ggplot2::element_blank(),
      plot.title   = ggplot2::element_text(hjust = 1)
    ) +
    ggplot2::scale_x_discrete(
      breaks = unique(region_data$net),
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25")
    ) +
    ggplot2::ggtitle(paste0("NIGHT - ", region))

  # Combine the day and night plots side by side
  combined_plot <- cowplot::plot_grid(
    day_plot, night_plot,
    align = "h",
    ncol  = 2
  )

  # Save the combined plot for the region
  ggplot2::ggsave(
    filename = paste0("combined_day_night_density_plot_", region, ".png"),
    plot     = combined_plot
  )

}

# Calculate density difference between day and night
density_difference <- Cal_data |>
  dplyr::filter(D.N %in% c("D", "N")) |>
  dplyr::group_by(Cluster, net) |>
  dplyr::summarise(
    day_density   = sum(density[D.N == "D"]),
    night_density = sum(density[D.N == "N"]),
    density_diff  = day_density - night_density
  ) |>
  dplyr::ungroup()

# Loop through each region to create Density Difference plots
for (region in unique(Cal_data$Region)) {
  
  # Filter the data by region
  region_data <- Cal_data |>
    dplyr::filter(Region == region)
  
  # Calculate density difference (Day - Night)
  density_difference <- region_data |>
    dplyr::group_by(net, Cluster) |>
    dplyr::summarize(
      day_density   = sum(density[D.N == "D"], na.rm = TRUE),
      night_density = sum(density[D.N == "N"], na.rm = TRUE)
    ) |>
    dplyr::mutate(density_diff = day_density - night_density)

  # Plot density difference
  density_diff_plot <- ggplot2::ggplot(
    data = density_difference,
    mapping = ggplot2::aes(
      x    = net,
      y    = density_diff,
      fill = as.factor(Cluster)
    )
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Depth (m)",
      y = "Density Difference (Day - Night) (ind m³)",
      fill = "Cluster"
    ) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_x_discrete(
      breaks = unique(density_difference$net),
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25")
    ) +
    ggplot2::ggtitle(paste0("Difference in Density (Day - Night) - ", region))

  # Save density difference plot for the region
  ggplot2::ggsave(
    filename = paste0("/scratch/srearl/density_difference_plot_", region, ".png"),
    plot     = density_diff_plot,
    width    = 10,
    height   = 6
  )

}

# Loop through each region to create Density per Cluster plots
for (region in unique(Cal_data$Region)) {
  
  # Filter the data by region
  region_data <- Cal_data |>
    dplyr::filter(Region == region)

  # Plot density per cluster
  dens_plot_cl <- ggplot2::ggplot(
    data = region_data,
    mapping = ggplot2::aes(
      x        = net,
      y        = density,
      group    = interaction(as.factor(Cluster), D.N),
      colour   = as.factor(Cluster),
      linetype = D.N,
      fill     = as.factor(Cluster),
      alpha    = D.N
    )
  ) +
    ggplot2::stat_summary(
      geom = "line",
      fun  = sum,
      size = 1.5
    ) +
    ggplot2::stat_summary(
      geom = "area",
      fun  = sum
    ) +
    ggplot2::facet_wrap(
      ~Cluster,
      scales = "free_y",
      nrow = 1
    ) +
    ggplot2::scale_x_discrete(
      limits = rev,
      breaks = unique(region_data$net),
      labels = c("1000", "800", "600", "400", "200", "100", "50", "25"),
      name   = "Depth (m)"
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_linetype_manual(
      values = c("solid", "dashed"),
      breaks = c("D", "N"),
      labels = c("Day", "Night"),
      name   = "Day/night"
    ) +
    ggplot2::scale_alpha_manual(
      values = c(0.1, 0.3),
      name   = "Day/night",
      breaks = c("D", "N"),
      labels = c("Day", "Night")
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      text = ggplot2::element_text(
        size = 25,
        face = "bold"
      ),
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill   = NA,
        size   = 3
      ),
      legend.position = "none",
      axis.text.y = ggplot2::element_text(
        angle = 90,
        hjust = 1
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1
      ),
      axis.title.x = ggplot2::element_text(
        angle = 180,
        vjust = 1
      ),
      plot.background = ggplot2::element_rect(fill = "transparent"),
      legend.key.width = ggplot2::unit(3, "line")
    ) +
    ggplot2::ylab(expression(bold("Density ind" ~ m^-3))) +
    ggplot2::guides(
      alpha    = ggplot2::guide_legend(nrow = 1),
      linetype = ggplot2::guide_legend(nrow = 1)
    )
  
  # Save density per cluster plot for the region
  ggplot2::ggsave(
    filename = paste0("/scratch/srearl/abund_plot_cal_", region, ".png"),
    plot     = dens_plot_cl,
    width    = 40,
    height   = 10,
    bg       = "white"
  )

}

####Plot transparency ~ size with marginal distribution#####
Cal_data$esd_mm  <- data_p$esd * .0053
Cal_data$Mean1   <- data_p$mean
Cal_data$Cluster <- as.factor(Cal_data$Cluster)

# Initialize a list to store the plots
plots_list <- list()

# Loop through each region
for (region in unique(Cal_data$Region)) {
  
  # Filter the data by region
  region_data <- Cal_data |>
    dplyr::filter(Region == region)
  # Summarize the data for the region
  df.mm <- region_data |>
    dplyr::group_by(Cluster) |>
    dplyr::summarise(
      MN.mean         = mean(Mean1),
      MN.esd          = mean(esd_mm),
      five_mean       = quantile(Mean1, probs = 0.05, na.rm = TRUE),
      nintyfive_mean  = quantile(Mean1, probs = 0.95, na.rm = TRUE),
      five_major      = quantile(esd_mm, probs = 0.05, na.rm = TRUE),
      nintyfive_major = quantile(esd_mm, probs = 0.95, na.rm = TRUE)
    )

  # Create hs plot for the region
  hs <- ggplot2::ggplot(
    data = df.mm,
    mapping = ggplot2::aes(
      x = MN.esd,
      y = MN.mean
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = Cluster),
      size = 5
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = five_mean,
        ymax = nintyfive_mean,
        color = Cluster
      ),
      width = 0.2,
      size = 2
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = five_major,
        xmax = nintyfive_major,
        color = Cluster
      ),
      height = 6,
      size = 2
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::ylab("Transparency") +
    ggplot2::xlab("Size class (mm)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 20),
      legend.position = "none",
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(
          t = 0,
          r = 21,
          b = 0,
          l = 0
        )
      )
    )
  
  # Create hs.mean1 plot for the region
  hs.mean1 <- ggplot2::ggplot(
    data = region_data,
    mapping = ggplot2::aes(
      x     = Mean1,
      y     = density,
      group = Cluster,
      color = Cluster,
      fill  = Cluster
    )
  ) +
    ggplot2::stat_summary_bin(
      geom  = "area",
      fun   = sum,
      bins  = 10,
      alpha = 0.3,
      size  = 1
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::theme_classic() +
    ggplot2::ylab(expression(bold("Density ind" ~ m^-3))) +
    ggpubr::clean_theme() +
    ggpubr::rremove("legend") +
    ggplot2::coord_flip()

  # Create hs.esd1 plot for the region
  hs.esd1 <- ggplot2::ggplot(
    data = region_data,
    mapping = ggplot2::aes(
      y     = density,
      x     = esd_mm,
      group = Cluster,
      color = Cluster,
      fill  = Cluster
    )
  ) +
    ggplot2::stat_summary_bin(
      geom  = "area",
      fun   = sum,
      bins  = 10,
      alpha = 0.3,
      size  = 1
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::theme_classic() +
    ggplot2::ylab(expression(bold("Density ind" ~ m^-3))) +
    ggpubr::clean_theme() +
    ggpubr::rremove("legend")
  
  # Insert hs.esd1 as the top x-axis and hs.mean1 as the right y-axis
  p1 <- cowplot::insert_xaxis_grob(hs, hs.esd1, grid::unit(.3, "null"), position = "top")
  p2 <- cowplot::insert_yaxis_grob(p1, hs.mean1, grid::unit(.3, "null"), position = "right")
  
  # Final combined plot for the region
  final_plot_cal <- cowplot::ggdraw(p2)
  
  # Save the plot in the list
  plots_list[[region]] <- final_plot_cal
  
  # Save the final plot for the region as an image
  ggplot2::ggsave(
    filename = paste0("/scratch/srearl/final_plot_cal_", region, ".png"),
    plot     = final_plot_cal,
    width    = 10,
    height   = 10,
    bg       = "white"
  )

}

# Display the list of plots
plots_list
