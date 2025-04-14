library(pointblank)

extract_columns <- function(
  df,
  debug = TRUE
  ) {
  
  pattern1 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"
  pattern2 <- "^([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([0-9]+)_([0-9]+_[0-9]+)$"
  pattern3 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+)$"
  pattern4 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[0-9]+)$"
  pattern5 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{4})_([0-9]+_[0-9]+)$"
  pattern6 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"
  pattern7 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_[a-zA-Z]_([0-9]+_[0-9]+)$"  

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
    yes  = base::sapply(matches1, function(x) x[2]),
    no   = df$cruise
  )
  df$moc <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[3]),
    no   = df$moc
  )
  df$net <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[4]),
    no   = df$net
  )
  df$fraction <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[5]),
    no   = df$fraction
  )
  df$lab_split <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[6]),
    no   = df$lab_split
  )
  df$split_fraction <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[7]),
    no   = df$split_fraction
  )
  df$photo_id <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[8]),
    no   = df$photo_id
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
  
  # extract components for pattern7
  matches7 <- base::regmatches(
    df$object_id,
    base::regexec(pattern7, df$object_id)
  )

  df$cruise <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes  = base::sapply(matches7, function(x) x[2]),
    no   = df$cruise
  )
  df$moc <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes  = base::sapply(matches7, function(x) x[3]),
    no   = df$moc
  )
  df$net <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes  = base::sapply(matches7, function(x) x[4]),
    no   = df$net
  )
  df$fraction <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes  = base::sapply(matches7, function(x) x[5]),
    no   = df$fraction
  )
  df$photo_id <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes  = base::sapply(matches7, function(x) x[6]),
    no   = df$photo_id
  )
  
  # perform pointblank test
  agent <- pointblank::create_agent(tbl = df) |>
    pointblank::col_vals_not_null(columns = dplyr::vars(cruise, moc, net, fraction)) |>
    pointblank::interrogate()

  if (debug == TRUE) {

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
    df$pattern <- base::ifelse(
      test = base::sapply(matches7, function(x) base::length(x) > 1),
      yes  = "pattern7",
      no   = df$pattern
    )

    cols <- c(
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
      
    df <- df |>
      dplyr::select(dplyr::any_of(cols))
    
  }
  
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


pattern4 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[0-9]+)$"
pattern6 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"

a problem is that both of these equate to pattern 6:
"209009_0200_1_512_1_1060"  <- matches both 4 and 6, needs to force it to be 4
"jc0214_m48_n9_2000_1_364"  <- matches6


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

base::regmatches(
  "10414_0000_01_1_20x_d_00080",
  base::regexec(pattern8, "10414_0000_01_1_20x_d_00080")
)
