## duck

CREATE TABLE agg AS 
SELECT * FROM read_csv(
  'eco_taxa.csv',
types={
  'process_time': 'VARCHAR',
  'acq_scan_time': 'VARCHAR',
  'acq_lut_16b_median': 'VARCHAR'
  },
nullstr = 'NA'
) ;

## -----

load_eco_taxa <- function(file_path) {

  ensure_numeric <- c(
    "object_area",
    "object_major",
    "object_minor",
    "object_esd",
    "object_depth_max",
    "object_depth_min",
    "acq_sub_part"
  )

  eco_taxa <- readr::read_delim(file_path)

  eco_taxa <- eco_taxa |>
    janitor::clean_names() |>
    tidyr::separate_wider_delim(
      col = object_id,
      delim = "_",
      names = c(
        "cruise",
        "moc",
        "net",
        "fraction"
      ),
      too_few = c("debug"),
      too_many = c("drop")
    ) |>
    dplyr::mutate(
      net = as.factor(net),
      dplyr::across(
        .cols = dplyr::all_of(ensure_numeric),
        .fns = as.numeric
      ),
      # do not convert cruise_moc_net to factor!
      cruise_moc_net = paste(
        cruise,
        moc,
        net,
        sep = "_"
      ),
      cruise_moc_net = tolower(cruise_moc_net),
      object_area_mm2 = dplyr::case_when(
        grepl(
          pattern = "4800",
          x = process_img_resolution,
          ignore.case = TRUE
        ) ~ object_area * (0.005291667^2),
        grepl(
          pattern = "2400",
          x = process_img_resolution,
          ignore.case = TRUE
        ) ~ object_area * (0.010583333^2),
        TRUE ~ NA_real_
      ),
      object_major_mm = dplyr::case_when(
        grepl(
          pattern = "4800",
          x = process_img_resolution,
          ignore.case = TRUE
        ) ~ object_major * 0.005291667,
        grepl(
          pattern = "2400",
          x = process_img_resolution,
          ignore.case = TRUE
        ) ~ object_major * 0.010583333,
        TRUE ~ NA_real_
      ),
      object_minor_mm = dplyr::case_when(
        grepl(
          pattern = "4800",
          x = process_img_resolution,
          ignore.case = TRUE
        ) ~ object_minor * 0.005291667,
        grepl(
          pattern = "2400",
          x = process_img_resolution,
          ignore.case = TRUE
        ) ~ object_minor * 0.010583333,
        TRUE ~ NA_real_
      ),
      object_esd_mm = dplyr::case_when(
        grepl(
          pattern = "4800",
          x = process_img_resolution,
          ignore.case = TRUE
        ) ~ object_esd * 0.005291667,
        grepl(
          pattern = "2400",
          x = process_img_resolution,
          ignore.case = TRUE
        ) ~ object_esd * 0.010583333,
        TRUE ~ NA_real_
      ),
      volume = (4 / 3) * pi * ((object_minor_mm * 0.5)^2) * (object_major_mm / 2),
      hdif = object_depth_max - object_depth_min, # SRE: are these the correct depths?
      split = (1 / acq_sub_part)
    )

  return(eco_taxa)
}

eco_taxa <- load_eco_taxa("~/Desktop/aggregates/ecotaxa_export_5446_20250307_1942.tsv")
eco_taxa <- load_eco_taxa("~/Desktop/gradients/ecotaxa_export_5421_20250307_2215.tsv")

moc_env <- readr::read_csv("~/localRepos/ciren/Amy_Gradients_MOCNESS_net_hydrography.csv") |>
  janitor::clean_names() |>
  dplyr::mutate(cruise_moc_net = tolower(cruise_moc_net))


# inner join will purge eco_taxa lacking env data
eco_env <- eco_taxa |>
  dplyr::inner_join(
    moc_env |>
      dplyr::select(
        -cruise,
        -net
      ),
    by = c("cruise_moc_net")
  )


eco_env <- eco_env |>
  dplyr::mutate(
    dry_weight = dplyr::case_when(
      grepl("Amphipoda", object_annotation_category, ignore.case = TRUE) ~ 0.0340 * volume,
      grepl("Calanoida", object_annotation_category, ignore.case = TRUE) ~ 0.0550 * volume,
      grepl("Chaetognatha", object_annotation_category, ignore.case = TRUE) ~ 0.0130 * volume,
      grepl("Decapoda", object_annotation_category, ignore.case = TRUE) ~ 0.0340 * volume,
      grepl("Euphausiacea", object_annotation_category, ignore.case = TRUE) ~ 0.0270 * volume,
      grepl("Foraminifera", object_annotation_category, ignore.case = TRUE) ~ 0.1420 * volume,
      grepl("Ostracoda", object_annotation_category, ignore.case = TRUE) ~ 0.0520 * volume,
      grepl("Poecilostomatoida", object_annotation_category, ignore.case = TRUE) ~ 0.0740 * volume,
      grepl("Thecosomata", object_annotation_category, ignore.case = TRUE) ~ 0.1913 * volume,
      TRUE ~ 0.0550 * volume # otherwise use Calanoida
      # TRUE ~ NA_real_
    ),
    o2_umol = (exp(-0.339 + (0.801 * log(dry_weight))) + 0.069 * (15)) / 22.4,
    co2 = (o2_umol) * 0.87 # o2 ~ co2 using a general rq
  )

eco_env <- eco_env |>
  # filter non-living
  dplyr::filter(
    !grepl(
      pattern = "not-living",
      x = object_annotation_hierarchy,
      ignore.case = TRUE
    )
  ) |>
  # add bins
  dplyr::mutate(
    bin = cut(
      x = object_esd_mm,
      breaks = c(
        seq(
          from = 0.25,
          to   = 74.25,
          by   = 0.25
        )
      ),
      labels = as.character(
        c(
          seq(
            from = 0.25,
            to   = 74,
            by   = 0.25
          )
        )
      )
    )
  )

summary_all <- eco_env |>
  dplyr::group_by(
    cruise,
    cruise_moc_net,
    station,
    d_n,
    net,
    bin
  ) |>
  dplyr::summarize(
    count           = dplyr::n(),
    depth_mean      = (mean(object_depth_min) + mean(object_depth_max)) / 2,
    depth_min       = min(object_depth_min),
    depth_max       = max(object_depth_max),
    hdif_median     = median(hdif), # SRE: each net should have the same value so is median needed?
    split_median    = median(split),
    frequency       = count / split_median, # use median split here and remainder?
    avg_sample_vol  = mean(sample_tot_vol), # from eco_taxa is this the correct volume?
    density_m3      = frequency / avg_sample_vol,
    # avg_bin_number  = mean(as.numeric(as.character(bin))),
    norm_bio_vol_m3 = (sum(volume) / split_median / avg_sample_vol),
    biomass_m3      = (sum(dry_weight) / split_median / avg_sample_vol),
    o2_m3           = (sum(o2_umol)) / avg_sample_vol / split_median,
    co2_m3          = (sum(co2) / avg_sample_vol / split_median),
    abundance_m2    = (frequency / avg_sample_vol * hdif_median),
    norm_bio_vol_m2 = norm_bio_vol_m3 * hdif_median,
    biomass_m2      = biomass_m3 * hdif_median,
    o2_m2           = o2_m3 * hdif_median,
    co2_m2          = co2_m3 * hdif_median
  ) |>
  dplyr::ungroup()


### BIOMASS SUMMARY
ggplot2::ggplot(
  data = summary_all,
  mapping = ggplot2::aes(
    x    = net,
    y    = biomass_m2,
    fill = "#D72000"
    # SRE: fraction was not included in the group by so how can it be used here?
    # alpha = fraction
  )
) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = c("#D72000")) +
  ggplot2::labs(
    x = "Net",
    y = expression("Biomass (Dry Weight)" ~ (mg ~ m^-2)),
    # SRE: the object H has not been defined in this workflow
    # title = paste(H, "Total Biomass by Net")
    title = "total biomass by net"
  ) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)) +
  ggplot2::scale_alpha_discrete(range = c(0.5, 1)) +
  ggplot2::guides(fill = FALSE) +
  ggplot2::ylim(0, 1000)

### OXYGEN USE SUMMARY
ggplot2::ggplot(
  data = summary_all,
  mapping = ggplot2::aes(
    x     = net,
    y     = o2_m2,
    fill  = "#FFAD0A",
    # SRE: fraction was not included in the grouping?
    # alpha = fraction
  )
) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = "#FFAD0A") +
  ggplot2::labs(
    x     = "Net",
    y     = expression(mu * mol ~ O[2] * m^-2 * h^-1),
    # SRE: the object H has not been defined in this workflow
    # title = paste(H, "Total Biomass by Net")
    title = "oxygen use by net"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face  = "bold",
      hjust = 0.5,
      size  = 11
    )
  ) +
  # SRE: scaling is extreme - definitely desired?
  ggplot2::scale_alpha_discrete(range = c(0.5, 1)) +
  ggplot2::guides(fill = FALSE) +
  ggplot2::ylim(0, 8000)

# WATERFALL
summary_all |>
  dplyr::group_by(
    net,
    bin
  ) |>
  dplyr::summarize(
    density_m3   = sum(density_m3),
    abundance_m2 = sum(abundance_m2)
  ) |>
  dplyr::mutate(bin_num = as.numeric(as.character(bin))) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x     = bin_num, # was bin but that is a factor, bin_num binN?
      y     = abundance_m2,
      color = net
    )
  ) +
  ggplot2::geom_point(size = 2) +
  # SRE: LaCroixColoR package is unavailable
  # ggplot2::scale_color_manual(
  #   values = (lacroix_palette("PeachPear", type = "continuous", n = 8)),
  #   labels = rev(net_labs)
  #   ) +
  ggplot2::labs(
    x     = expression("Size class" ~ (mm^3)),
    y     = expression("Abundance" ~ (particles ~ m^-2)),
    title = paste("Oct 2018 (Day)"),
    color = ""
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(reverse = T)) +
  ggplot2::scale_y_log10(
    limits = c(0.001, 100000), # you may need to change your scale
    breaks = c(0.1, 1, 10, 100, 1000, 10000, 100000), # you may need to change your scale
    labels = c("0.1", "1", "10", "100", "1000", "10000", "100000")
  ) +
  # you can limit your scale based on what biomass you effectively sample, but
  # you should look to see all data first
  ggplot2::scale_x_log10(
    limits = c(.001, 1000),
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
    labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000")
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face  = "bold",
      hjust = 0.5,
      size  = 14
    ),
    legend.position = "bottom"
  ) +
  ggplot2::theme(legend.text = ggplot2::element_text(size = 11)) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    axis.text.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 14),
    strip.text = ggplot2::element_text(size = 14, face = "bold")
  )


summary_all |>
  dplyr::mutate(bin_num = as.numeric(as.character(bin))) |>
  dplyr::filter(bin_num > 0.01 & bin_num < 100) |>
  dplyr::group_by(net) |>
  dplyr::group_split(.keep = TRUE) |>
  purrr::map(
    .f = ~ summary(
      lm(
        formula = log(abundance_m2) ~ log(bin_num),
        data    = .x,
        y       = TRUE
      )
    )
  ) |>
  capture.output(file = "/tmp/stats_abundance_bin.txt")

# DAY NIGHT

day_night_sept_2022 <- summary_all |>
  dplyr::mutate(bin_num = as.numeric(as.character(bin))) |>
  dplyr::filter(
    grepl(
      pattern     = "nh1208",
      x           = cruise_moc_net,
      ignore.case = TRUE
    ),
    grepl(
      pattern     = "_m7|_m8",
      x           = cruise_moc_net,
      ignore.case = TRUE
    ),
    bin_num >= 0.01 & bin_num <= 100
  ) |>
  dplyr::group_by(
    d_n,
    net
  ) |>
  dplyr::summarise(
    depth_mean      = sum(depth_mean),
    density_m3      = sum(density_m3),
    norm_bio_vol_m3 = sum(norm_bio_vol_m3),
    biomass_m3      = sum(biomass_m3),
    o2_m3           = sum(o2_m3),
    co2_m3          = sum(co2_m3),
    abundance_m2    = sum(abundance_m2),
    norm_bio_vol_m2 = sum(norm_bio_vol_m2),
    biomass_m2      = sum(biomass_m2),
    o2_m2           = sum(o2_m2),
    co2_m2          = sum(co2_m2),
  ) |>
  dplyr::ungroup()

str(day_night_sept_2022)

numeric_columns <- day_night_sept_2022 |>
  dplyr::select(
    tidyselect::where(is.numeric)
  ) |>
  colnames()

estimate_day_night <- function(variable) {
  variable_sym <- rlang::sym(variable)

  if (grepl("depth", variable_sym, ignore.case = TRUE)) {
    day_night_calcs <- day_night_sept_2022 |>
      dplyr::select(
        d_n,
        net,
        !!variable_sym
      ) |>
      tidyr::pivot_wider(
        names_from  = d_n,
        values_from = !!variable_sym
      ) |>
      dplyr::mutate(
        depth_mean = (D + N) / 2
      ) |>
      # below for LONG
      dplyr::select(
        -D,
        -N
      ) |>
      tidyr::expand_grid(
        MR = c("migratory", "resident")
      )
    # below for WIDE
    # dplyr::select(
    #   -D,
    #   -N
    # )
  } else {
    day_night_calcs <- day_night_sept_2022 |>
      dplyr::select(
        d_n,
        net,
        !!variable_sym
      ) |>
      tidyr::pivot_wider(
        names_from  = d_n,
        values_from = !!variable_sym
      ) |>
      dplyr::mutate(
        migratory = abs(D - N),
        resident  = pmin(D, N)
      ) |>
      # below for LONG
      dplyr::select(
        -D,
        -N
      ) |>
      tidyr::pivot_longer(
        cols = c(
          migratory,
          resident
        ),
        names_to = "MR",
        values_to = variable
      )
    # below for WIDE
    # dplyr::rename(
    #   !!paste0(variable, "_mig") := migratory,
    #   !!paste0(variable, "_res") := resident
    # ) |>
    # dplyr::select(
    #   -D,
    #   -N
    # )
  }

  return(day_night_calcs)
}

# just one
# estimate_day_night("norm_bio_vol_m2")

# Example usage with a list of column names
# variables <- c("norm_bio_vol_m2", "biomass_m2", "o2_m2", "co2_m2")

day_night_calcs <- purrr::map(
  .x = numeric_columns,
  .f = estimate_day_night
) |>
  purrr::reduce(
    .f = dplyr::left_join,
    by = c("net", "MR")
  )

# just one plot mirroring Amy's plot
day_night_calcs |>
  ggplot2::ggplot(
    ggplot2::aes(
      x     = as.factor(depth_mean),
      y     = norm_bio_vol_m2,
      fill  = "grey7",
      alpha = MR
    )
  ) +
  ggplot2::scale_fill_manual(values = ("grey7")) +
  ggplot2::scale_alpha_discrete(
    range  = c(0.5, 1),
    labels = c("migratory", "resident"),
    drop   = FALSE
  ) +
  ggplot2::geom_col(
    position = "stack",
    na.rm = FALSE
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x     = "minimum net depth (m)",
    y     = expression("mg Biomass" ~ m^-2),
    title = "Dry Weight Biomass",
    alpha = ""
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face  = "bold",
      hjust = 0.5,
      size  = 12
    )
  ) +
  ggplot2::guides(fill = FALSE) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line        = ggplot2::element_line(colour = "black"),
    axis.text.x      = ggplot2::element_text(size = 12),
    axis.text.y      = ggplot2::element_text(size = 12, vjust = -1),
    axis.ticks.y     = (ggplot2::element_blank()),
    axis.title       = ggplot2::element_text(size = 14),
    strip.text       = ggplot2::element_text(size = 14, face = "bold"),
    legend.text      = ggplot2::element_text(size = 11),
    legend.position  = "bottom"
  )


# Define the numerical variables to be plotted
numerical_vars <- c(
  "density_m3",
  "norm_bio_vol_m3",
  "biomass_m3",
  "o2_m3",
  "co2_m3",
  "abundance_m2",
  "norm_bio_vol_m2",
  "biomass_m2",
  "o2_m2",
  "co2_m2"
)

# Convert dataset to long format for faceting
day_night_long <- day_night_calcs |>
  tidyr::pivot_longer(
    cols = tidyselect::all_of(numerical_vars),
    names_to = "variable",
    values_to = "value"
  )


# Plot all variables in a faceted grid
ggplot2::ggplot(
  data = day_night_long,
  mapping = ggplot2::aes(
    x    = value,
    y    = as.factor(depth_mean),
    fill = MR
  )
) +
  ggplot2::geom_col(
    position = "stack",
    na.rm = FALSE
  ) + # Stacked bars
  ggplot2::scale_fill_manual(values = c("migratory" = "gray40", "resident" = "gray70")) + # Gray shades
  ggplot2::facet_wrap(~variable, scales = "free_x") + # Faceted grid, free x-axis scale
  ggplot2::labs(
    x     = "value",
    y     = "minimum net depth (m)",
    title = "stacked bar plots of environmental variables",
    fill  = "movement type"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.line        = ggplot2::element_line(colour = "black"),
    axis.text.x      = ggplot2::element_text(size = 12),
    axis.text.y      = ggplot2::element_text(size = 12),
    axis.ticks.y     = ggplot2::element_blank(),
    axis.title       = ggplot2::element_text(size = 14),
    strip.text       = ggplot2::element_text(size = 14, face = "bold"),
    legend.text      = ggplot2::element_text(size = 11),
    legend.position  = "bottom"
  )


# HEAT MAPS -------------

summary_all |>
  dplyr::mutate(bin_num = as.numeric(as.character(bin))) |>
  dplyr::filter(
    grepl(
      pattern     = "nh1208",
      x           = cruise_moc_net,
      ignore.case = TRUE
    ),
    grepl(
      pattern     = "_m7|_m8",
      x           = cruise_moc_net,
      ignore.case = TRUE
    )
    # bin_num >= 0.01 & bin_num <= 100
  ) |>
  dplyr::group_by(
    bin,
    net
  ) |>
  dplyr::summarise(
    depth_mean      = sum(depth_mean),
    density_m3      = sum(density_m3),
    norm_bio_vol_m3 = sum(norm_bio_vol_m3),
    biomass_m3      = sum(biomass_m3),
    o2_m3           = sum(o2_m3),
    co2_m3          = sum(co2_m3),
    abundance_m2    = sum(abundance_m2),
    norm_bio_vol_m2 = sum(norm_bio_vol_m2),
    biomass_m2      = sum(biomass_m2),
    o2_m2           = sum(o2_m2),
    co2_m2          = sum(co2_m2),
  ) |>
  dplyr::ungroup()


# SCRATCH ---------------

day <- day_night_sept_2022 |>
  dplyr::filter(d_n == "D") |>
  dplyr::rename(
    Tot_BV_m2 = norm_bio_vol_m2,
    Tot_Ox_m2 = o2_m2,
    med_depth = depth_mean
  ) |>
  dplyr::select(
    d_n,
    net,
    Tot_BV_m2,
    Tot_Ox_m2,
    med_depth
  )

night <- day_night_sept_2022 |>
  dplyr::filter(d_n == "N") |>
  dplyr::rename(
    Tot_BV_m2 = norm_bio_vol_m2,
    Tot_Ox_m2 = o2_m2,
    med_depth = depth_mean
  ) |>
  dplyr::select(
    d_n,
    net,
    Tot_BV_m2,
    Tot_Ox_m2,
    med_depth
  )

# This is the standard code, but use it always
DayNight <- data.frame(
  # "BV_Mig" = c(abs(day$Tot_BV_m2 - night$Tot_BV_m2)),
  "BV_Mig" = abs(day$Tot_BV_m2 - night$Tot_BV_m2),
  # "BV_Res" = do.call(pmin, (as.data.frame(cbind(day$Tot_BV_m2, night$Tot_BV_m2)))),
  "BV_Res" = pmin(day$Tot_BV_m2, night$Tot_BV_m2),
  # "BM_Mig" = c(abs(day$Tot_BM_m2 - night$Tot_BM_m2)),
  # "BM_Res" = do.call(pmin, (as.data.frame(cbind(day$Tot_BM_m2, night$Tot_BM_m2)))),
  "Ox_Mig" = c(abs(day$Tot_Ox_m2 - night$Tot_Ox_m2)),
  "Ox_Res" = do.call(pmin, (as.data.frame(cbind(day$Tot_Ox_m2, night$Tot_Ox_m2)))),
  # "CO2_Mig" = c(abs(day$Tot_CO2_m2 - night$Tot_CO2_m2)),
  # "CO2_Mig" = do.call(pmin, (as.data.frame(cbind(day$Tot_CO2_m2, night$Tot_CO2_m2)))),
  # "BM_DVM" = c((day$Tot_BM_m2 - night$Tot_BM_m2)),
  # "BM_day" = c(day$Tot_BM_m2),
  # "BM_night" = c(night$Tot_BM_m2),
  # "BV_DVM" = c((day$Tot_BV_m2 - night$Tot_BV_m2)),
  "Med_Depth" = as.factor(apply(as.data.frame(cbind(day$med_depth, night$med_depth)), 1, FUN = mean)),
  "Net" = c(day$net)
)


DN2 <- data.frame(
  "Net" = as.factor(c(1:8, 1:8)),
  "M_R" = as.factor(c(rep.int("M", 8), rep.int("R", 8))),
  # "BM" = c(DayNight$BM_Mig, DayNight$BM_Res),
  "BV" = c(DayNight$BV_Mig, DayNight$BV_Res),
  "Ox" = c(DayNight$Ox_Mig, DayNight$Ox_Res),
  # "CO2" = c(DayNight$CO2_Mig, DayNight$CO2_Res),
  "Med_Depth" = as.factor(rep.int(DayNight$Med_Depth, 2))
)

identical(alpha, beta)
identical(beta, charlie)
