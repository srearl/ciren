library(pointblank)

extract_columns <- function(
  ecotaxa_file,
  debug = TRUE
  ) {
  
  # initialize new columns with NA
  ecotaxa_file$cruise         <- NA
  ecotaxa_file$photo_id       <- NA

flowcam_pattern <- "^([0-9]{5})_([0-9]{4})_([0-9]{2})_([0-9]{1})_([0-9]+[a-zA-Z]+)_([a-zA-Z])_([0-9]+)$"

# randomly select a subsample of 100 records
set.seed(123) # set seed for reproducibility
subsample <- ecotaxa_file[sample(nrow(ecotaxa_file), 100), ]

# check if all rows in the subsample match the flowcam_pattern
matches_flowcam_pattern <- base::regmatches(
  subsample$object_id,
  base::regexec(flowcam_pattern, subsample$object_id)
)

# test if all rows match the pattern
flowcam_pattern_true <- all(sapply(matches_flowcam_pattern, function(x) length(x) > 1))

# message(flowcam_pattern_true)
# flowcam_pattern_true <- TRUE

if (flowcam_pattern_true == TRUE) {

ecotaxa_file <- extract_flowcam_columns(
  ecotaxa_file = ecotaxa_file,
  pattern      = flowcam_pattern,
  debug        = FALSE
)

} else {

test_cols <- c(
  "cruise",
  "photo_id",
  "moc",
  "net",
  "fraction",
  "lab_split",
  "split_fraction"
)

  ecotaxa_file$moc            <- NA
  ecotaxa_file$net            <- NA
  ecotaxa_file$fraction       <- NA
  ecotaxa_file$lab_split      <- NA
  ecotaxa_file$split_fraction <- NA
  
  pattern1 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"
  pattern2 <- "^([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([0-9]+)_([0-9]+_[0-9]+)$"
  pattern3 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+)$"
  pattern4 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]+)_([0-9]+_[0-9]+)_([0-9]+_[0-9]+)$"
  pattern5 <- "^([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{4})_([0-9]+_[0-9]+)$"
  pattern6 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z0-9]+)_([0-9]+_[0-9]+)$"
  pattern7 <- "^([a-zA-Z]+[0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_([a-zA-Z][0-9]+)_[a-zA-Z]_([0-9]+_[0-9]+)$"  
  
  # extract components for pattern1
  matches1 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern1, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[2]),
    no   = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[3]),
    no   = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[4]),
    no   = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[5]),
    no   = ecotaxa_file$fraction
  )
  ecotaxa_file$lab_split <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[6]),
    no   = ecotaxa_file$lab_split
  )
  ecotaxa_file$split_fraction <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[7]),
    no   = ecotaxa_file$split_fraction
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches1, function(x) base::length(x) > 1),
    yes  = base::sapply(matches1, function(x) x[8]),
    no   = ecotaxa_file$photo_id
  )
  
  # extract components for pattern2
  matches2 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern2, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches2, function(x) base::length(x) > 1),
    yes  = base::sapply(matches2, function(x) x[2]),
    no   = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches2, function(x) base::length(x) > 1),
    yes  = base::sapply(matches2, function(x) x[3]),
    no   = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches2, function(x) base::length(x) > 1),
    yes  = base::sapply(matches2, function(x) x[4]),
    no   = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches2, function(x) base::length(x) > 1),
    yes  = base::sapply(matches2, function(x) x[5]),
    no   = ecotaxa_file$fraction
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches2, function(x) base::length(x) > 1),
    yes  = base::sapply(matches2, function(x) x[6]),
    no   = ecotaxa_file$photo_id
  )
  
  # extract components for pattern3
  matches3 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern3, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes  = base::sapply(matches3, function(x) x[2]),
    no   = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes  = base::sapply(matches3, function(x) x[3]),
    no   = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes  = base::sapply(matches3, function(x) x[4]),
    no   = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes  = base::sapply(matches3, function(x) x[5]),
    no   = ecotaxa_file$fraction
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches3, function(x) base::length(x) > 1),
    yes  = base::sapply(matches3, function(x) x[6]),
    no   = ecotaxa_file$photo_id
  )
  
  # extract components for pattern4
  matches4 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern4, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[2]),
    no   = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[3]),
    no   = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[4]),
    no   = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[5]),
    no   = ecotaxa_file$fraction
  )
  ecotaxa_file$lab_split <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[6]),
    no   = ecotaxa_file$lab_split
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches4, function(x) base::length(x) > 1),
    yes  = base::sapply(matches4, function(x) x[7]),
    no   = ecotaxa_file$photo_id
  )
  
  # extract components for pattern5
  matches5 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern5, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes  = base::sapply(matches5, function(x) x[2]),
    no   = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes  = base::sapply(matches5, function(x) x[3]),
    no   = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes  = base::sapply(matches5, function(x) x[4]),
    no   = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes  = base::sapply(matches5, function(x) x[5]),
    no   = ecotaxa_file$fraction
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches5, function(x) base::length(x) > 1),
    yes  = base::sapply(matches5, function(x) x[6]),
    no   = ecotaxa_file$photo_id
  )

  # extract components for pattern6
  matches6 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern6, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes  = base::sapply(matches6, function(x) x[2]),
    no   = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes  = base::sapply(matches6, function(x) x[3]),
    no   = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes  = base::sapply(matches6, function(x) x[4]),
    no   = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes  = base::sapply(matches6, function(x) x[5]),
    no   = ecotaxa_file$fraction
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches6, function(x) base::length(x) > 1),
    yes  = base::sapply(matches6, function(x) x[6]),
    no   = ecotaxa_file$photo_id
  )
  
  # extract components for pattern7
  matches7 <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern7, ecotaxa_file$object_id)
  )

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes  = base::sapply(matches7, function(x) x[2]),
    no   = ecotaxa_file$cruise
  )
  ecotaxa_file$moc <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes  = base::sapply(matches7, function(x) x[3]),
    no   = ecotaxa_file$moc
  )
  ecotaxa_file$net <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes  = base::sapply(matches7, function(x) x[4]),
    no   = ecotaxa_file$net
  )
  ecotaxa_file$fraction <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes  = base::sapply(matches7, function(x) x[5]),
    no   = ecotaxa_file$fraction
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches7, function(x) base::length(x) > 1),
    yes  = base::sapply(matches7, function(x) x[6]),
    no   = ecotaxa_file$photo_id
  )
  
}
  
  # perform pointblank test
  # agent <- pointblank::create_agent(tbl = ecotaxa_file) |>
  #   pointblank::col_vals_not_null(columns = dplyr::vars(cruise, moc, net, fraction)) |>
  #   pointblank::interrogate()

  if (debug == TRUE) {

    print(agent)
    
    ecotaxa_file$pattern <- NA

    ecotaxa_file$pattern <- base::ifelse(
      test = base::sapply(matches1, function(x) base::length(x) > 1),
      yes  = "pattern1",
      no   = ecotaxa_file$pattern
    )
    ecotaxa_file$pattern <- base::ifelse(
      test = base::sapply(matches2, function(x) base::length(x) > 1),
      yes  = "pattern2",
      no   = ecotaxa_file$pattern
    )
    ecotaxa_file$pattern <- base::ifelse(
      test = base::sapply(matches3, function(x) base::length(x) > 1),
      yes  = "pattern3",
      no   = ecotaxa_file$pattern
    )
    ecotaxa_file$pattern <- base::ifelse(
      test = base::sapply(matches4, function(x) base::length(x) > 1),
      yes  = "pattern4",
      no   = ecotaxa_file$pattern
    )
    ecotaxa_file$pattern <- base::ifelse(
      test = base::sapply(matches5, function(x) base::length(x) > 1),
      yes  = "pattern5",
      no   = ecotaxa_file$pattern
    )
    ecotaxa_file$pattern <- base::ifelse(
      test = base::sapply(matches6, function(x) base::length(x) > 1),
      yes  = "pattern6",
      no   = ecotaxa_file$pattern
    )
    ecotaxa_file$pattern <- base::ifelse(
      test = base::sapply(matches7, function(x) base::length(x) > 1),
      yes  = "pattern7",
      no   = ecotaxa_file$pattern
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
      
    ecotaxa_file <- ecotaxa_file |>
      dplyr::select(dplyr::any_of(cols))
    
  }
  
  return(ecotaxa_file)

}


extract_flowcam_columns <- function(
  ecotaxa_file,
  pattern = flowcam_pattern,
  debug = FALSE
) {
  
  # temporary
  # ecotaxa_file$cruise         <- NA
  # ecotaxa_file$photo_id       <- NA
  
  test_cols <- c(
    "cruise",
    "photo_id",
    "depth",
    "niskin",
    "mode",
    "magnification",
    "duplicates_removed"
  )

  # extract components for pattern1
  matches_flowcam <- base::regmatches(
    ecotaxa_file$object_id,
    base::regexec(pattern, ecotaxa_file$object_id)
  )

  ecotaxa_file$depth              <- NA
  ecotaxa_file$niskin             <- NA
  ecotaxa_file$mode               <- NA
  ecotaxa_file$magnification      <- NA
  ecotaxa_file$duplicates_removed <- NA

  ecotaxa_file$cruise <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes  = base::sapply(matches_flowcam, function(x) x[2]),
    no   = ecotaxa_file$cruise
  )
  ecotaxa_file$depth <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes  = base::sapply(matches_flowcam, function(x) x[3]),
    no   = ecotaxa_file$depth
  )
  ecotaxa_file$niskin <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes  = base::sapply(matches_flowcam, function(x) x[4]),
    no   = ecotaxa_file$niskin
  )
  ecotaxa_file$mode <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes  = base::sapply(matches_flowcam, function(x) x[5]),
    no   = ecotaxa_file$mode
  )
  ecotaxa_file$magnification <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes  = base::sapply(matches_flowcam, function(x) x[6]),
    no   = ecotaxa_file$magnification
  )
  ecotaxa_file$duplicates_removed <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes  = base::sapply(matches_flowcam, function(x) x[7]),
    no   = ecotaxa_file$duplicates_removed
  )
  ecotaxa_file$photo_id <- base::ifelse(
    test = base::sapply(matches_flowcam, function(x) base::length(x) > 1),
    yes  = base::sapply(matches_flowcam, function(x) x[8]),
    no   = ecotaxa_file$photo_id
  )

return(ecotaxa_file)

}

# TESTING -----

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

# load flowcam
loaded <- read.delim(
    file             = "~/Desktop/dataset_var/Dataset_variations/Rhizaria_Flowcam_14986_20250205_1901.tsv",
    sep              = "\t",
    header           = TRUE,
    stringsAsFactors = FALSE
    )

extracted <- extract_columns(loaded, debug = FALSE)

flowcam_pattern <- "^([0-9]{5})_([0-9]{4})_([0-9]{2})_([0-9]{1})_([0-9]+[a-zA-Z]+)_([a-zA-Z])_([0-9]+)$"
extracted <- extract_flowcam_columns(loaded, pattern = flowcam_pattern, debug = FALSE)

"10414_0000_01_1_20x_d_00080"


# SCRATCH ------

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
    flowcam_pattern, c(
      "10414_0000_01_1_20x_d_00080",
      "10414_0000_01_1_20x_d_00090",
    "not"
      )
)
)
