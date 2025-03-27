library(pointblank)

extract_columns <- function(df) {

  # define the regex patterns
  # pattern1 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([a-zA-Z0-9]+)_([0-9]+)_([0-9]+)$"
  # pattern2 <- "^([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)$"
  # pattern3 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+)$"
  # pattern4 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)_([0-9]+)$"
  # pattern5 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{4})_([0-9]+)_([0-9]+)$"
  
  pattern1 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"
  pattern2 <- "^([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([0-9]+)_([0-9]+_[0-9]+)$"
  pattern3 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+)$"
  pattern4 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[0-9]+)$"
  pattern5 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{4})_([0-9]+_[0-9]+)$"
  pattern6 <- "^([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"
  
  # initialize new columns with NA
  df$cruise         <- NA
  df$moc            <- NA
  df$net            <- NA
  df$fraction       <- NA
  df$lab_split      <- NA
  df$split_fraction <- NA
  df$photo_id       <- NA
  
  # extract components for pattern1
  matches1 <- base::regmatches(
    df$object_id,
    base::regexec(pattern1, df$object_id)
  )

  df$cruise <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[2]),
    no = df$cruise
  )
  df$moc <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[3]),
    no = df$moc
  )
  df$net <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[4]),
    no = df$net
  )
  df$fraction <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[5]),
    no = df$fraction
  )
  df$lab_split <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[6]),
    no = df$lab_split
  )
  df$split_fraction <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[7]),
    no = df$split_fraction
  )
  df$photo_id <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes = base::sapply(matches1, function(x) x[8]),
    no = df$photo_id
  )
  
  # extract components for pattern2
  matches2 <- base::regmatches(
    df$object_id,
    base::regexec(pattern2, df$object_id)
  )

  df$cruise <- base::ifelse(
    test = base::sapply(matches2, function(x) base::length(x) > 1),
    yes  = base::sapply(matches2, function(x) x[2]),
    no   = df$cruise
  )
  df$moc <- base::ifelse(
    test = base::sapply(matches2, function(x) base::length(x) > 1),
    yes  = base::sapply(matches2, function(x) x[3]),
    no   = df$moc
  )
  df$net <- base::ifelse(
    test = base::sapply(matches2, function(x) base::length(x) > 1),
    yes  = base::sapply(matches2, function(x) x[4]),
    no   = df$net
  )
  df$fraction <- base::ifelse(
    test = base::sapply(matches2, function(x) base::length(x) > 1),
    yes  = base::sapply(matches2, function(x) x[5]),
    no   = df$fraction
  )
  df$photo_id <- base::ifelse(
    test = base::sapply(matches2, function(x) base::length(x) > 1),
    yes  = base::sapply(matches2, function(x) x[6]),
    no   = df$photo_id
  )
  
  # extract components for pattern3
  matches3 <- base::regmatches(
    df$object_id,
    base::regexec(pattern3, df$object_id)
  )

  df$cruise <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes  = base::sapply(matches3, function(x) x[2]),
    no   = df$cruise
  )
  df$moc <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes  = base::sapply(matches3, function(x) x[3]),
    no   = df$moc
  )
  df$net <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes  = base::sapply(matches3, function(x) x[4]),
    no   = df$net
  )
  df$fraction <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes  = base::sapply(matches3, function(x) x[5]),
    no   = df$fraction
  )
  df$photo_id <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes  = base::sapply(matches3, function(x) x[6]),
    no   = df$photo_id
  )
  
  # extract components for pattern4
  matches4 <- base::regmatches(
    df$object_id,
    base::regexec(pattern4, df$object_id)
  )

  df$cruise <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[2]),
    no   = df$cruise
  )
  df$moc <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[3]),
    no   = df$moc
  )
  df$net <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[4]),
    no   = df$net
  )
  df$fraction <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[5]),
    no   = df$fraction
  )
  df$lab_split <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[6]),
    no   = df$lab_split
  )
  df$photo_id <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[7]),
    no   = df$photo_id
  )
  
  # extract components for pattern5
  matches5 <- base::regmatches(
    df$object_id,
    base::regexec(pattern5, df$object_id)
  )

  df$cruise <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes  = base::sapply(matches5, function(x) x[2]),
    no   = df$cruise
  )
  df$moc <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes  = base::sapply(matches5, function(x) x[3]),
    no   = df$moc
  )
  df$net <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes  = base::sapply(matches5, function(x) x[4]),
    no   = df$net
  )
  df$fraction <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes  = base::sapply(matches5, function(x) x[5]),
    no   = df$fraction
  )
  df$photo_id <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes  = base::sapply(matches5, function(x) x[6]),
    no   = df$photo_id
  )

  # extract components for pattern6
  matches6 <- base::regmatches(
    df$object_id,
    base::regexec(pattern6, df$object_id)
  )

  df$cruise <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes  = base::sapply(matches6, function(x) x[2]),
    no   = df$cruise
  )
  df$moc <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes  = base::sapply(matches6, function(x) x[3]),
    no   = df$moc
  )
  df$net <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes  = base::sapply(matches6, function(x) x[4]),
    no   = df$net
  )
  df$fraction <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes  = base::sapply(matches6, function(x) x[5]),
    no   = df$fraction
  )
  df$photo_id <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes  = base::sapply(matches6, function(x) x[6]),
    no   = df$photo_id
  )
  
  # Perform pointblank test
  agent <- pointblank::create_agent(tbl = df) |>
    pointblank::col_vals_not_null(columns = dplyr::vars(cruise, moc, net, fraction)) |>
    pointblank::interrogate()

  # FOR DEVELOPMENT ONLY -----

  print(agent)
  
  df$pattern <- NA

  df$pattern <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = "pattern1",
    no   = df$pattern
  )
  df$pattern <- base::ifelse(
    test = base::sapply(matches2, function(x) base::length(x) > 1),
    yes  = "pattern2",
    no   = df$pattern
  )
  df$pattern <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes  = "pattern3",
    no   = df$pattern
  )
  df$pattern <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = "pattern4",
    no   = df$pattern
  )
  df$pattern <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes  = "pattern5",
    no   = df$pattern
  )
  df$pattern <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes  = "pattern6",
    no   = df$pattern
  )

  cols <- c("object_id", "cruise", "moc", "net", "fraction", "lab_split", "split_fraction", "photo_id", "pattern")
    
  df <- df |>
    dplyr::select(dplyr::any_of(cols))
  
  # END DEVELOPMENT ONLY -----
  
  return(df)

}

# Load the data
loaded <- read.delim(
    file             = "~/Desktop/dataset_var/Dataset_variations/ecotaxa_export_NA2021_MOCNESS.tsv",
    sep              = "\t",
    header           = TRUE,
    stringsAsFactors = FALSE
    )

# Apply the extract_columns function
extracted <- extract_columns(loaded)

# Load another dataset and apply the function
loaded <- read.delim(
    file             = "~/Desktop/aggregates/ecotaxa_export_5446_20250307_1942.tsv",
    sep              = "\t",
    header           = TRUE,
    stringsAsFactors = FALSE
    )

extracted <- extract_columns(loaded)

# a problem with pattern6:
# ae2112_m22_n2_d1_1_75
# ae2112_m22_n2_d2_a_1_1


# Load another dataset and apply the function
loaded <- read.delim(
    file             = "~/Desktop/gradients/ecotaxa_export_5421_20250307_2215_rm.tsv",
    sep              = "\t",
    header           = TRUE,
    stringsAsFactors = FALSE
    )

extracted <- extract_columns(loaded)

# same problem with pattern6 as above

# test_object_ids <- c(
#   "204707_0200_1_1",
#   "204707_0200_1_32_1_6",
#   "204707_0200_1_128_1_tot_1_1",
#   "201009_0200_1_256_2_tot_1_1"
# )

# matches <- base::regmatches(test_object_ids, base::regexec(pattern1, test_object_ids))
# print(matches)
