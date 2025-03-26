#' @title Load and minimally process downloaded EcoTaxa data
#'
#' @description This function loads and processes EcoTaxa data from a TSV file.
#'   It ensures that the data is properly formatted and annotated for further
#'   analysis.
#'
#' @param file_path (character) Quoted path to the EcoTaxa TSV file.
#'
#' @return A data frame with the processed EcoTaxa data, including additional
#'   columns for cruise, moc, net, fraction, and select object dimensions
#'   converted from pixels to millimeters.
#'
#' @importFrom tools file_ext
#' @importFrom readr read_delim
#' @importFrom janitor clean_names
#' @importFrom tidyr separate_wider_delim
#' @importFrom dplyr mutate across all_of case_when group_by summarize filter
#' @importFrom pointblank create_agent col_vals_not_null col_vals_lte
#'   col_vals_equal interrogate
#'
#' @export
#'
#' @examples \dontrun{ eco_taxa_df <- load_eco_taxa("path/to/eco_taxa_file.tsv")
#'   print(eco_taxa_df) }
#'
load_eco_taxa <- function(file_path) {

  if (tools::file_ext(file_path) != "tsv") {
    stop("provided file must be of type tsv")
  }

  eco_taxa <- read.delim(
    file             = file_path,
    sep              = "\t",
    header           = TRUE,
    stringsAsFactors = FALSE
    )
  # eco_taxa <- readr::read_delim(file_path)

  ensure_numeric <- c(
    "object_area",
    "object_major",
    "object_minor",
    "object_esd",
    "object_depth_max",
    "object_depth_min",
    "acq_sub_part"
  )

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
      too_few  = c("debug"),
      too_many = c("drop")
    ) |>
    dplyr::mutate(
      # revisit converting net to factor
      # net = as.factor(net),
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
        grepl("4800", process_img_resolution, ignore.case = TRUE) ~ object_area * (0.005291667^2),
        grepl("2400", process_img_resolution, ignore.case = TRUE) ~ object_area * (0.010583333^2),
        TRUE ~ NA_real_
      ),
      object_major_mm = dplyr::case_when(
        grepl("4800", process_img_resolution, ignore.case = TRUE) ~ object_major * 0.005291667,
        grepl("2400", process_img_resolution, ignore.case = TRUE) ~ object_major * 0.010583333,
        TRUE ~ NA_real_
      ),
      object_minor_mm = dplyr::case_when(
        grepl("4800", process_img_resolution, ignore.case = TRUE) ~ object_minor * 0.005291667,
        grepl("2400", process_img_resolution, ignore.case = TRUE) ~ object_minor * 0.010583333,
        TRUE ~ NA_real_
      ),
      object_esd_mm = dplyr::case_when(
        grepl("4800", process_img_resolution, ignore.case = TRUE) ~ object_esd * 0.005291667,
        grepl("2400", process_img_resolution, ignore.case = TRUE) ~ object_esd * 0.010583333,
        TRUE ~ NA_real_
      ),
      volume = (4 / 3) * pi * ((object_minor_mm * 0.5)^2) * (object_major_mm / 2),
      hdif = object_depth_max - object_depth_min,
      split = (1 / acq_sub_part)
    )

  # data validations:
  # 1. check for null cruise_moc_net
  # 2. check for multiple dates per cruise-moc pair
  # 3. check for inconsistent latitudes per cruise-moc pair
  agent <- pointblank::create_agent(tbl = eco_taxa) |>
    pointblank::col_vals_not_null(columns = vars(cruise_moc_net)) |>
    pointblank::col_vals_lte(
      columns = vars(unique_dates),
      value = 1,
      preconditions = function(x) {
        x |>
          dplyr::group_by(
            cruise,
            moc
          ) |>
          dplyr::summarize(
            unique_dates = dplyr::n_distinct(object_date),
            .groups      = "drop"
          )
      }
    ) |>
    pointblank::col_vals_equal(
      columns = vars(stdev),
      value = 0,
      preconditions = function(x) {
        x |>
          dplyr::group_by(
            cruise,
            moc
          ) |>
          dplyr::summarize(
            stdev   = sd(object_lat, na.rm = TRUE),
            .groups = "drop"
          )
      }
    ) |>
    pointblank::interrogate()


  if (agent$validation_set$all_passed[1] == FALSE) {

    cruise_moc_net_null <- eco_taxa |>
      dplyr::filter(is.na(cruise_moc_net))

    warning(
      "encountered NULL cruise_moc_net(s): ",
      cruise_moc_net_null
    )

  } 
  
  if (agent$validation_set$all_passed[2] == FALSE) {

    cruise_moc_date <- eco_taxa |>
      dplyr::group_by(
        cruise,
        moc
      ) |>
      dplyr::summarize(
        unique_dates = dplyr::n_distinct(object_date),
        .groups      = "drop"
        ) |>
      dplyr::filter(unique_dates != 1)

    stop(
      "encountered a cruise-moc pair with more than one date: ",
      paste(
        cruise_moc_date$cruise,
        cruise_moc_date$moc,
        collapse = " "
      )
    )

  }

  if (agent$validation_set$all_passed[3] == FALSE) {

    cruise_moc_lat <- eco_taxa |>
      dplyr::group_by(
        cruise,
        moc
      ) |>
      dplyr::summarize(
        stdev   = sd(object_lat, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::filter(stdev != 0)

    warning(
      "encountered a cruise-moc pair with inconsistent latitude(s): ",
      paste(
        cruise_moc_lat$cruise,
        cruise_moc_lat$moc,
        collapse = " "
      )
    )

  }

  eco_taxa <- annotate_daytime_parallel(eco_taxa)

  return(eco_taxa)

}