

extract_columns <- function(
  ecotaxa_file,
  debug = FALSE
) {
  # randomly select a subsample of 1000 records
  set.seed(123)
  subsample <- ecotaxa_file[sample(nrow(ecotaxa_file), 1000), ]

  flowcam_pattern_true <- FALSE
  moc_pattern_true <- FALSE
  uvp_pattern_true <- FALSE

  # FLOWCAM

  flowcam_pattern <- "^([0-9]{5})_([0-9]{4})_([0-9]{2})_([0-9]{1})_([0-9]+[a-zA-Z]+)_([a-zA-Z])_([0-9]+)$"
  # flowcam_pattern: "10414_0000_01_1_20x_d_00080"

  flowcam_pattern2 <- "^([0-9]{5})_([0-9]{4})_([0-9]{2})_([0-9]{1})_([0-9]+[a-zA-Z]+)_([0-9]+_[a-zA-Z])_([0-9]+)$"
  # flowcam_pattern2:  "10423_0800_22_1_20x_2_d_00116"

  # check if rows in subsample match flowcam_patterns
  matches_flowcam_pattern <- base::regmatches(
    subsample$object_id,
    base::regexec(flowcam_pattern, subsample$object_id)
  )
  matches_flowcam_pattern2 <- base::regmatches(
    subsample$object_id,
    base::regexec(flowcam_pattern2, subsample$object_id)
  )

  flowcam_pattern_true <- all(
    sapply(matches_flowcam_pattern, function(x) length(x) > 1) |
      sapply(matches_flowcam_pattern2, function(x) length(x) > 1)
  )

  if (flowcam_pattern_true == TRUE) {
    message("flowcam pattern matched")

    flowcam_patterns <- list(
      flowcam_pattern = flowcam_pattern,
      flowcam_pattern2 = flowcam_pattern2
    )

    extracted_file <- extract_flowcam_columns(
      ecotaxa_file = ecotaxa_file,
      pattern = flowcam_patterns,
      debug = debug
    )

    return(extracted_file)
  }

  # MOC

  pattern1 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{4})_([0-9]+_[0-9]+)_([a-zA-Z0-9]+_[a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"
  pattern3 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+)$"
  pattern4 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{4})_([0-9]{1}_[0-9]+)_([0-9]{1}_[0-9]+)$"
  pattern5 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{4})_([0-9]+_[0-9]+)$"
  pattern6 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"
  pattern7 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_[a-zA-Z]_([0-9]+_[0-9]+)$"
  # pattern2 <- "^([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([0-9]+)_([0-9]+_[0-9]+)$"

  matches1 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern1, ecotaxa_file$object_id)
  )
  matches3 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern3, ecotaxa_file$object_id)
  )
  matches4 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern4, ecotaxa_file$object_id)
  )
  matches5 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern5, ecotaxa_file$object_id)
  )
  matches6 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern6, ecotaxa_file$object_id)
  )
  matches7 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern7, ecotaxa_file$object_id)
  )

  moc_pattern_true <- all(
    sapply(matches1, function(x) length(x) > 1) |
      sapply(matches3, function(x) length(x) > 1) |
      sapply(matches4, function(x) length(x) > 1) |
      sapply(matches5, function(x) length(x) > 1) |
      sapply(matches6, function(x) length(x) > 1) |
      sapply(matches7, function(x) length(x) > 1)
  )

  if (moc_pattern_true == TRUE) {
    message("MOC pattern matched")

    moc_patterns <- list(
      pattern1 = pattern1,
      pattern3 = pattern3,
      pattern4 = pattern4,
      pattern5 = pattern5,
      pattern6 = pattern6,
      pattern7 = pattern7
    )

    extracted_file <- extract_moc_columns(
      ecotaxa_file = ecotaxa_file,
      pattern = moc_patterns,
      debug = debug
    )

    return(extracted_file)
  }

  # UVP

  uvp_pattern <- "^([0-9]{4})([0-9]{2})([0-9]{2})-([0-9]{2})([0-9]{2})([0-9]{2})-([0-9]{3})_([0-9]+)$"

  # check if rows in subsample match uvp_pattern
  matches_uvp_pattern <- base::regmatches(
    subsample$object_id,
    base::regexec(uvp_pattern, subsample$object_id)
  )

  uvp_pattern_true <- all(
    sapply(matches_uvp_pattern, function(x) length(x) > 1)
  )

  if (uvp_pattern_true == TRUE) {
    message("uvp pattern matched")

    uvp_patterns <- list(
      uvp_patterns = uvp_pattern
    )

    extracted_file <- extract_uvp_columns(
      ecotaxa_file = ecotaxa_file,
      pattern = uvp_patterns,
      debug = debug
    )

    return(extracted_file)
  }

  if (
    !flowcam_pattern_true &&
      !moc_pattern_true &&
      !uvp_pattern_true
  ) {
    message("could match the object_id to any patterns")
    return(NULL)
  }
}

extract_moc_columns <- function(
  ecotaxa_file,
  pattern = moc_patterns,
  debug = FALSE
) {
  ecotaxa_file$cruise <- NA
  ecotaxa_file$photo_id <- NA
  ecotaxa_file$moc <- NA
  ecotaxa_file$net <- NA
  ecotaxa_file$fraction <- NA
  ecotaxa_file$lab_split <- NA
  ecotaxa_file$split_fraction <- NA

  # extract components for pattern1
  matches1 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern$pattern1, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[2]),
    no = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[3]),
    no = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[4]),
    no = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[5]),
    no = ecotaxa_file$fraction
  )
  ecotaxa_file$lab_split <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[6]),
    no = ecotaxa_file$lab_split
  )
  ecotaxa_file$split_fraction <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[7]),
    no = ecotaxa_file$split_fraction
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[8]),
    no = ecotaxa_file$photo_id
  )

  # extract components for pattern3
  matches3 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern$pattern3, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes = base::sapply(matches3, function(x) x[2]),
    no = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes = base::sapply(matches3, function(x) x[3]),
    no = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes = base::sapply(matches3, function(x) x[4]),
    no = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes = base::sapply(matches3, function(x) x[5]),
    no = ecotaxa_file$fraction
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes = base::sapply(matches3, function(x) x[6]),
    no = ecotaxa_file$photo_id
  )

  # extract components for pattern4
  matches4 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern$pattern4, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes = base::sapply(matches4, function(x) x[2]),
    no = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes = base::sapply(matches4, function(x) x[3]),
    no = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes = base::sapply(matches4, function(x) x[4]),
    no = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes = base::sapply(matches4, function(x) x[5]),
    no = ecotaxa_file$fraction
  )
  ecotaxa_file$lab_split <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes = base::sapply(matches4, function(x) x[6]),
    no = ecotaxa_file$lab_split
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes = base::sapply(matches4, function(x) x[7]),
    no = ecotaxa_file$photo_id
  )

  # extract components for pattern5
  matches5 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern$pattern5, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes = base::sapply(matches5, function(x) x[2]),
    no = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes = base::sapply(matches5, function(x) x[3]),
    no = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes = base::sapply(matches5, function(x) x[4]),
    no = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes = base::sapply(matches5, function(x) x[5]),
    no = ecotaxa_file$fraction
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes = base::sapply(matches5, function(x) x[6]),
    no = ecotaxa_file$photo_id
  )

  # extract components for pattern6
  matches6 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern$pattern6, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes = base::sapply(matches6, function(x) x[2]),
    no = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes = base::sapply(matches6, function(x) x[3]),
    no = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes = base::sapply(matches6, function(x) x[4]),
    no = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes = base::sapply(matches6, function(x) x[5]),
    no = ecotaxa_file$fraction
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes = base::sapply(matches6, function(x) x[6]),
    no = ecotaxa_file$photo_id
  )

  # extract components for pattern7
  matches7 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern$pattern7, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes = base::sapply(matches7, function(x) x[2]),
    no = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes = base::sapply(matches7, function(x) x[3]),
    no = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes = base::sapply(matches7, function(x) x[4]),
    no = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes = base::sapply(matches7, function(x) x[5]),
    no = ecotaxa_file$fraction
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes = base::sapply(matches7, function(x) x[6]),
    no = ecotaxa_file$photo_id
  )

  # add the pattern invoked for testing and debugging (optional)

  patterns <- list(
    pattern$pattern1,
    pattern$pattern3,
    pattern$pattern4,
    pattern$pattern5,
    pattern$pattern6,
    pattern$pattern7
  )

  pattern_names <- c(
    "pattern1",
    "pattern3",
    "pattern4",
    "pattern5",
    "pattern6",
    "pattern7"
  )

  # use purrr::map2() to iterate over patterns and pattern_names
  ecotaxa_file$pattern <- NA
  ecotaxa_file$pattern <- purrr::reduce(
    .x = purrr::map2(
      patterns,
      pattern_names,
      ~ {
        base::ifelse(
          test = base::sapply(
            base::regmatches(
              ecotaxa_file$object_id,
              base::regexec(.x, ecotaxa_file$object_id)
            ),
            function(x) base::length(x) > 1
          ),
          yes = .y,
          no = ecotaxa_file$pattern
        )
      }
    ),
    .f = ~ ifelse(
      test = is.na(.x),
      yes = .y,
      no = .x
    )
  )

  agent <- pointblank::create_agent(tbl = ecotaxa_file) |>
    pointblank::col_vals_not_null(
      columns = dplyr::vars(
        cruise,
        moc,
        net,
        fraction
      ),
      label = "cruise_moc_net_fraction"
    ) |>
    pointblank::col_vals_not_null(
      columns = c("lab_split"),
      preconditions = function(x) {
        x |>
          dplyr::filter(pattern %in% c("pattern1", "pattern4"))
      }
    ) |>
    pointblank::col_vals_not_null(
      columns = c("split_fraction"),
      preconditions = function(x) {
        x |>
          dplyr::filter(pattern %in% c("pattern1"))
      }
    ) |>
    pointblank::interrogate()

  if (any(agent$validation_set$all_passed == FALSE)) {
    failed_indices <- which(agent$validation_set$all_passed == FALSE)

    purrr::walk(
      .x = failed_indices,
      .f = function(i) {
        failed_brief <- agent$validation_set$brief[i]
        failed_rows <- agent$validation_set$n_failed[i]
        message(glue::glue(
          "validation failed for: {failed_brief} ({failed_rows} rows failed)"
        ))
      }
    )
  }

  moc_cols <- c(
    "object_id",
    "cruise",
    "moc",
    "net",
    "fraction",
    "lab_split",
    "split_fraction",
    "photo_id",
    "pattern"
  )

  if (debug == TRUE) {
    print(agent)

    ecotaxa_file <- ecotaxa_file |>
      dplyr::select(dplyr::any_of(moc_cols))
  } else {
    # remove any parsed columns that are empty
    ecotaxa_file <- ecotaxa_file |>
      dplyr::select(
        -dplyr::any_of(moc_cols[sapply(
          ecotaxa_file[moc_cols],
          function(x) all(is.na(x))
        )])
      ) |>
      dplyr::select(-pattern) # remove pattern column
  }

  return(ecotaxa_file)
}


extract_flowcam_columns <- function(
  ecotaxa_file,
  pattern = flowcam_patterns,
  debug = FALSE
) {
  ecotaxa_file$cruise <- NA
  ecotaxa_file$photo_id <- NA
  ecotaxa_file$depth <- NA
  ecotaxa_file$niskin <- NA
  ecotaxa_file$mode <- NA
  ecotaxa_file$magnification <- NA
  ecotaxa_file$duplicates_removed <- NA

  # extract components for flowcam_pattern
  matches_flowcam <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern$flowcam_pattern, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam, function(x) x[2]),
    no = ecotaxa_file$cruise
  )
  ecotaxa_file$depth <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam, function(x) x[3]),
    no = ecotaxa_file$depth
  )
  ecotaxa_file$niskin <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam, function(x) x[4]),
    no = ecotaxa_file$niskin
  )
  ecotaxa_file$mode <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam, function(x) x[5]),
    no = ecotaxa_file$mode
  )
  ecotaxa_file$magnification <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam, function(x) x[6]),
    no = ecotaxa_file$magnification
  )
  ecotaxa_file$duplicates_removed <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam, function(x) x[7]),
    no = ecotaxa_file$duplicates_removed
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam, function(x) x[8]),
    no = ecotaxa_file$photo_id
  )

  # extract components for flowcam_pattern2
  matches_flowcam2 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern$flowcam_pattern2, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches_flowcam2, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam2, function(x) x[2]),
    no = ecotaxa_file$cruise
  )
  ecotaxa_file$depth <- base::ifelse(
    test = base::sapply(matches_flowcam2, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam2, function(x) x[3]),
    no = ecotaxa_file$depth
  )
  ecotaxa_file$niskin <- base::ifelse(
    test = base::sapply(matches_flowcam2, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam2, function(x) x[4]),
    no = ecotaxa_file$niskin
  )
  ecotaxa_file$mode <- base::ifelse(
    test = base::sapply(matches_flowcam2, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam2, function(x) x[5]),
    no = ecotaxa_file$mode
  )
  ecotaxa_file$magnification <- base::ifelse(
    test = base::sapply(matches_flowcam2, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam2, function(x) x[6]),
    no = ecotaxa_file$magnification
  )
  ecotaxa_file$duplicates_removed <- base::ifelse(
    test = base::sapply(matches_flowcam2, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam2, function(x) x[7]),
    no = ecotaxa_file$duplicates_removed
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches_flowcam2, function(x) base::length(x) > 1),
    yes = base::sapply(matches_flowcam2, function(x) x[8]),
    no = ecotaxa_file$photo_id
  )

  patterns <- list(
    pattern$flowcam_pattern,
    pattern$flowcam_pattern2
  )

  pattern_names <- c(
    "flowcam_pattern",
    "flowcam_pattern2"
  )

  # use purrr::map2() to iterate over patterns and pattern_names
  ecotaxa_file$pattern <- NA
  ecotaxa_file$pattern <- purrr::reduce(
    .x = purrr::map2(
      patterns,
      pattern_names,
      ~ {
        base::ifelse(
          test = base::sapply(
            base::regmatches(
              ecotaxa_file$object_id,
              base::regexec(.x, ecotaxa_file$object_id)
            ),
            function(x) base::length(x) > 1
          ),
          yes = .y,
          no = ecotaxa_file$pattern
        )
      }
    ),
    .f = ~ ifelse(
      test = is.na(.x),
      yes = .y,
      no = .x
    )
  )

  agent <- pointblank::create_agent(tbl = ecotaxa_file) |>
    pointblank::col_vals_not_null(
      columns = dplyr::vars(
        cruise,
        photo_id,
        depth,
        niskin,
        mode,
        magnification,
        duplicates_removed
      ),
      label = "flowcam_cols"
    ) |>
    pointblank::interrogate()

  if (any(agent$validation_set$all_passed == FALSE)) {
    failed_indices <- which(agent$validation_set$all_passed == FALSE)

    purrr::walk(
      .x = failed_indices,
      .f = function(i) {
        failed_brief <- agent$validation_set$brief[i]
        failed_rows <- agent$validation_set$n_failed[i]
        message(glue::glue(
          "validation failed for: {failed_brief} ({failed_rows} rows failed)"
        ))
      }
    )
  }

  flowcam_cols <- c(
    "object_id",
    "cruise",
    "depth",
    "niskin",
    "mode",
    "magnification",
    "duplicates_removed",
    "photo_id",
    "pattern"
  )

  if (debug == TRUE) {
    print(agent)

    ecotaxa_file <- ecotaxa_file |>
      dplyr::select(dplyr::any_of(flowcam_cols))
  } else {
    ecotaxa_file <- ecotaxa_file |>
      # remove any parsed columns that are empty
      dplyr::select(
        -dplyr::any_of(flowcam_cols[sapply(
          ecotaxa_file[flowcam_cols],
          function(x) all(is.na(x))
        )])
      ) |>
      # remove pattern column
      dplyr::select(-pattern)
  }

  return(ecotaxa_file)
}

extract_uvp_columns <- function(
  ecotaxa_file,
  pattern = uvp_patterns,
  debug = FALSE
) {
  ecotaxa_file$year <- NA
  ecotaxa_file$month <- NA
  ecotaxa_file$day <- NA
  ecotaxa_file$hour <- NA
  ecotaxa_file$minute <- NA
  ecotaxa_file$second <- NA
  ecotaxa_file$millisecond <- NA
  ecotaxa_file$object_number <- NA

  # extract components for flowcam_pattern
  matches_uvp <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern$uvp_pattern, ecotaxa_file$object_id)
  )

  ecotaxa_file$year <- base::ifelse(
    test = base::sapply(matches_uvp, function(x) base::length(x) > 1),
    yes = base::sapply(matches_uvp, function(x) x[2]),
    no = ecotaxa_file$year
  )
  ecotaxa_file$month <- base::ifelse(
    test = base::sapply(matches_uvp, function(x) base::length(x) > 1),
    yes = base::sapply(matches_uvp, function(x) x[3]),
    no = ecotaxa_file$month
  )
  ecotaxa_file$day <- base::ifelse(
    test = base::sapply(matches_uvp, function(x) base::length(x) > 1),
    yes = base::sapply(matches_uvp, function(x) x[4]),
    no = ecotaxa_file$day
  )
  ecotaxa_file$hour <- base::ifelse(
    test = base::sapply(matches_uvp, function(x) base::length(x) > 1),
    yes = base::sapply(matches_uvp, function(x) x[5]),
    no = ecotaxa_file$hour
  )
  ecotaxa_file$minute <- base::ifelse(
    test = base::sapply(matches_uvp, function(x) base::length(x) > 1),
    yes = base::sapply(matches_uvp, function(x) x[6]),
    no = ecotaxa_file$minute
  )
  ecotaxa_file$second <- base::ifelse(
    test = base::sapply(matches_uvp, function(x) base::length(x) > 1),
    yes = base::sapply(matches_uvp, function(x) x[7]),
    no = ecotaxa_file$second
  )
  ecotaxa_file$millisecond <- base::ifelse(
    test = base::sapply(matches_uvp, function(x) base::length(x) > 1),
    yes = base::sapply(matches_uvp, function(x) x[8]),
    no = ecotaxa_file$millisecond
  )
  ecotaxa_file$object_number <- base::ifelse(
    test = base::sapply(matches_uvp, function(x) base::length(x) > 1),
    yes = base::sapply(matches_uvp, function(x) x[9]),
    no = ecotaxa_file$object_number
  )

  patterns <- list(
    pattern$uvp_pattern
  )

  pattern_names <- c(
    "uvp_pattern"
  )

  # use purrr::map2() to iterate over patterns and pattern_names
  ecotaxa_file$pattern <- NA
  ecotaxa_file$pattern <- purrr::reduce(
    .x = purrr::map2(
      patterns,
      pattern_names,
      ~ {
        base::ifelse(
          test = base::sapply(
            base::regmatches(
              ecotaxa_file$object_id,
              base::regexec(.x, ecotaxa_file$object_id)
            ),
            function(x) base::length(x) > 1
          ),
          yes = .y,
          no = ecotaxa_file$pattern
        )
      }
    ),
    .f = ~ ifelse(
      test = is.na(.x),
      yes = .y,
      no = .x
    )
  )

  agent <- pointblank::create_agent(tbl = ecotaxa_file) |>
    pointblank::col_vals_not_null(
      columns = dplyr::vars(
        year,
        month,
        day,
        hour,
        minute,
        second,
        millisecond,
        object_number
      ),
      label = "uvp_cols"
    ) |>
    pointblank::interrogate()

  if (any(agent$validation_set$all_passed == FALSE)) {
    failed_indices <- which(agent$validation_set$all_passed == FALSE)

    purrr::walk(
      .x = failed_indices,
      .f = function(i) {
        failed_brief <- agent$validation_set$brief[i]
        failed_rows <- agent$validation_set$n_failed[i]
        message(glue::glue(
          "validation failed for: {failed_brief} ({failed_rows} rows failed)"
        ))
      }
    )
  }

  uvp_cols <- c(
    "object_id",
    "year",
    "month",
    "day",
    "hour",
    "minute",
    "second",
    "millisecond",
    "object_number"
  )

  if (debug == TRUE) {
    print(agent)

    ecotaxa_file <- ecotaxa_file |>
      dplyr::select(dplyr::any_of(uvp_cols))
  } else {
    ecotaxa_file <- ecotaxa_file |>
      # remove any parsed columns that are empty
      dplyr::select(
        -dplyr::any_of(uvp_cols[sapply(
          ecotaxa_file[uvp_cols],
          function(x) all(is.na(x))
        )])
      ) |>
      # remove pattern column
      dplyr::select(-pattern)
  }

  return(ecotaxa_file)
}

# TESTING -----

# test: ecotaxa_export_NP2018_MOCNESS.tsv

pattern1 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{4})_([0-9]+_[0-9]+)_([a-zA-Z0-9]+_[a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"

# cmn    frac lbsp spfr  foto_id
"101005_5000_1_2_1a_tot_1_1"
"101002_0200_1_128_1_tot_1_1"

base::regmatches(
  c(
    "101005_5000_1_2_1a_tot_1_1",
    "101002_0200_1_128_1_tot_1_1",
    "not"
  ),
  base::regexec(
    pattern1,
    c(
      "101005_5000_1_2_1a_tot_1_1",
      "101002_0200_1_128_1_tot_1_1",
      "not"
    )
  )
)

loaded <- read.delim(
  file = "~/Desktop/dataset_var/Dataset_variations/ecotaxa_export_NP2018_MOCNESS.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

extracted <- extract_columns(loaded, debug = TRUE)


# test: ecotaxa_export_NA2021_MOCNESS.tsv

loaded <- read.delim(
  file = "~/Desktop/dataset_var/Dataset_variations/ecotaxa_export_NA2021_MOCNESS.tsv",
  # file             = "/tmp/ecotaxa_export_NA2021_MOCNESS.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

extracted <- extract_columns(loaded, debug = TRUE)


# test: ecotaxa_export_5446_20250307_1942.tsv

loaded <- read.delim(
  file = "~/Desktop/aggregates/ecotaxa_export_5446_20250307_1942.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

extracted <- extract_columns(loaded, debug = FALSE)

# a problem with pattern6:
# ae2112_m22_n2_d1_1_75
# ae2112_m22_n2_d2_a_1_1

# test: ecotaxa_export_5421_20250307_2215_rm.tsv
loaded <- read.delim(
  file = "~/Desktop/gradients/ecotaxa_export_5421_20250307_2215_rm.tsv",
  # file = "~/Desktop/gradients/ecotaxa_export_5421_20250307_2215.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

extracted <- extract_columns(loaded, debug = FALSE)

# same problem with pattern6 as above

# test_object_ids <- c(
#   "204707_0200_1_1",
#   "204707_0200_1_32_1_6",
#   "204707_0200_1_128_1_tot_1_1",
#   "201009_0200_1_256_2_tot_1_1"
# )

# matches <- base::regmatches(test_object_ids, base::regexec(pattern1, test_object_ids))
# print(matches)

# load flowcam
loaded <- read.delim(
  file = "~/Desktop/dataset_var/Dataset_variations/Rhizaria_Flowcam_14986_20250205_1901.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

extracted <- extract_columns(loaded, debug = FALSE)
readr::write_csv(extracted[c(1:1000, 4000:5000), ], "/tmp/flow.csv")

# extracted <- extract_flowcam_columns(
#   loaded,
#   pattern = flowcam_pattern,
#   debug = FALSE
# )

flowcam_pattern <- "^([0-9]{5})_([0-9]{4})_([0-9]{2})_([0-9]{1})_([0-9]+[a-zA-Z]+)_([a-zA-Z])_([0-9]+)$"
flowcam_pattern2 <- "^([0-9]{5})_([0-9]{4})_([0-9]{2})_([0-9]{1})_([0-9]+[a-zA-Z]+)_([0-9]+_[a-zA-Z])_([0-9]+)$"

"10414_0000_01_1_20x_d_00080"
"10423_0800_22_1_20x_2_d_00116"

base::regmatches(
  c(
    "10414_0000_01_1_20x_d_00080",
    "10423_0800_22_1_20x_2_d_00116",
    "not"
  ),
  base::regexec(
    pattern7,
    c(
      "10414_0000_01_1_20x_d_00080",
      "10423_0800_22_1_20x_2_d_00116",
      "not"
    )
  )
)

# UVP

"20230629-193335-678_1"

uvp_pattern <- "^([0-9]{4})([0-9]{2})([0-9]{2})-([0-9]{2})([0-9]{2})([0-9]{2})-([0-9]{3})_([0-9]+)$"

base::regmatches(
  c(
    "20230629-193335-678_1",
    "not"
  ),
  base::regexec(
    uvp_pattern,
    c(
      "20230629-193335-678_1",
      "not"
    )
  )
)

loaded <- read.delim(
  file = "~/Desktop/dataset_var/Dataset_variations/BATS_UVP_13454_20250204_1932.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

extracted <- extract_columns(loaded, debug = TRUE)
extracted <- extract_columns(loaded, debug = FALSE)


# SCRATCH ------

pattern4 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[0-9]+)$"
pattern6 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"

# a problem is that both of these equate to pattern 6:
"209009_0200_1_512_1_1060" # matches both 4 and 6, needs to force it to be 4
"jc0214_m48_n9_2000_1_364" # matches6


base::regmatches(
  "209009_0200_1_512_1_1060",
  base::regexec(pattern4, "209009_0200_1_512_1_1060")
)

base::regmatches(
  "209009_0200_1_512_1_1060",
  base::regexec(pattern6, "209009_0200_1_512_1_1060")
)

base::regmatches(
  "jc0214_m48_n9_2000_1_364",
  base::regexec(pattern4, "jc0214_m48_n9_2000_1_364")
)

base::regmatches(
  "jc0214_m48_n9_2000_1_364",
  base::regexec(pattern6, "jc0214_m48_n9_2000_1_364")
)

#                 ae2306            m35               n8               d3               b     1_1235
pattern7 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_[a-zA-Z]_([0-9]+_[0-9]+)$"

"ae2112_m22_n1_d1_1_14"
"ae2306_m35_n8_d3_b_1_1235"

base::regmatches(
  "ae2306_m35_n8_d3_b_1_1235",
  base::regexec(pattern6, "ae2306_m35_n8_d3_b_1_1235")
)

base::regmatches(
  "ae2112_m22_n1_d1_1_14",
  base::regexec(pattern6, "ae2112_m22_n1_d1_1_14")
)

base::regmatches(
  "ae2306_m35_n8_d3_b_1_1235",
  base::regexec(pattern7, "ae2306_m35_n8_d3_b_1_1235")
)

base::regmatches(
  "ae2112_m22_n1_d1_1_14",
  base::regexec(pattern7, "ae2112_m22_n1_d1_1_14")
)


# FLOWCAM

# FlowCam is:
# 5-digit BATS cruise number,
# 4-digit depth,
# two digit niskin number,
# 1 for AutoImage mode or 2 for trigger mode,
# magnification,
# d (meaning duplicate particles were removed),
# and last number is the image number!

pattern8 <- "^([0-9]{5})_([0-9]{4})_([0-9]{2})_([0-9]{1})_([0-9]+[a-zA-Z]+)_([a-zA-Z])_([0-9]+)$"

"10414_0000_01_1_20x_d_00080"
"10414_0000_01_1_20x_d_00090"

base::regmatches(
  "10414_0000_01_1_20x_d_00080",
  base::regexec(flowcam_pattern, "10414_0000_01_1_20x_d_00080")
)

base::regmatches(
  c(
    "10414_0000_01_1_20x_d_00080",
    "10414_0000_01_1_20x_d_00090",
    "not"
  ),
  base::regexec(
    flowcam_pattern,
    c(
      "10414_0000_01_1_20x_d_00080",
      "10414_0000_01_1_20x_d_00090",
      "not"
    )
  )
)

# pattern 4 problem

# one
# 201009_0200_1_256_2_tot_1_1

# four
# 209010_0200_1_1024_1_1068

# cmn    frac lb_sp  foto_id
# 209010_0200_1_1024_1_1068
# 2,3,4  5    6      7

pattern4 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[0-9]+)$"

pattern1 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"
pattern2 <- "^([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([0-9]+)_([0-9]+_[0-9]+)$"
pattern3 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+)$"
pattern4 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[0-9]+)$"
pattern5 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{4})_([0-9]+_[0-9]+)$"
pattern6 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"
pattern7 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_[a-zA-Z]_([0-9]+_[0-9]+)$"

base::regmatches(
  c(
    "201009_0200_1_256_2_tot_1_1",
    "209010_0200_1_1024_1_1068",
    "not"
  ),
  base::regexec(
    pattern4,
    c(
      "201009_0200_1_256_2_tot_1_1",
      "209010_0200_1_1024_1_1068",
      "not"
    )
  )
)


# pattern 2 versus 4

# solution was to omit 2

# cruise moc net frac foto_id
"jc0214_m47_n10_0200_1_4"

# cmn     frac lb_sp foto_id
"204702_0200_1_32_1_1"

base::regmatches(
  c(
    "jc0214_m47_n10_0200_1_4", # this is pattern 2 and 6
    "204702_0200_1_32_1_1", # this is pattern 2 and 4
    "jc0214_m48_n9_5000_1_1319", # this is pattern 2 and 6
    "not"
  ),
  base::regexec(
    pattern6,
    c(
      "jc0214_m47_n10_0200_1_4",
      "204702_0200_1_32_1_1",
      "jc0214_m48_n9_5000_1_1319",
      "not"
    )
  )
)

pattern2 <- "^([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([0-9]+)_([0-9]+_[0-9]+)$"
pattern4 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{4})_([0-9]{1}_[0-9]+)_([0-9]{1}_[0-9]+)$"
pattern6 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"

pattern1 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"
pattern3 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+)$"
pattern5 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{4})_([0-9]+_[0-9]+)$"
pattern7 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_[a-zA-Z]_([0-9]+_[0-9]+)$"

# pattern 4 versus 6

"209009_0200_1_512_1_1060"

base::regmatches(
  c(
    "jc0214_m47_n10_0200_1_4", # this is pattern 2 and 6
    "204702_0200_1_32_1_1", # this is pattern 2 and 4
    "jc0214_m48_n9_5000_1_1319", # this is pattern 2 and 6
    "209009_0200_1_512_1_1060", # this is pattern 2 and 4
    "not"
  ),
  base::regexec(
    pattern7,
    c(
      "jc0214_m47_n10_0200_1_4",
      "204702_0200_1_32_1_1",
      "jc0214_m48_n9_5000_1_1319",
      "209009_0200_1_512_1_1060",
      "not"
    )
  )
)
