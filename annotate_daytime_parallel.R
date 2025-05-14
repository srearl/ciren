#' @title Annotate EcoTaxa records as day or night based on ship position, and
#'   collection date and time.
#'
#' @description \code{annotate_daytime_parallel} annotates whether a given
#'   observation in an EcoTaxa record was collected during the day or night
#'   based on ship position, and collection date and time. The function assumes
#'   that the collection time is provided in UTC. Minimal processing of the
#'   EcoTaxa file such that cruise and MOC ids are included as separate columns
#'   in the data frame as `cruise` and `moc`, respectively, is required.
#'
#' @note The function was parallelized as the initial implementation computed
#'   `is_day` for each record of the data frame. Since all MOCs do (or should)
#'   have the same collection details, the function was refactored to group the
#'   records by cruise and MOC, and computes `is_day` for each group, which are
#'   then joined to the full data set. This reduces the number of calls to
#'   `is_daytime` and speeds up the computation, but is still computationally
#'   expensive and, thus, the parallelization was retained.
#'
#' @param eco_taxa_df (character) Unquoted name of data downloaded from EcoTaxa
#'   to the R environment at least minimally processed to have separate columns
#'   for the cruise and MOC ids as `cruise` and `moc`, respectively.
#' @param workers (integer) The number of parallel workers to use. Defaults to
#'   one less than the available cores.
#'
#' @return A data frame with an additional column `is_day` indicating whether
#'   the observation occurred during daytime.
#'
#' @importFrom dplyr group_by summarize left_join select ungroup mutate n
#' @importFrom future availableCores plan multisession
#' @importFrom furrr future_pmap_lgl
#' @importFrom SunCalcMeeus is_daytime
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' eco_taxa_df <- data.frame(
#'   cruise = c("cruise1", "cruise1"),
#'   moc = c("moc1", "moc1"),
#'   object_lat = c(34.5, 34.5),
#'   object_lon = c(-120.5, -120.5),
#'   object_date = c("2025-03-26", "2025-03-26"),
#'   object_time = c("12:00:00", "13:00:00")
#' )
#' result <- annotate_daytime_parallel(eco_taxa_df)
#' print(result)
#' }
#'
#' @export
#'

eco_taxa_df <- this_file

annotate_daytime_parallel <- function(
  eco_taxa_df,
  workers = future::availableCores() - 1
) {
  required_cols <- c(
    "object_lat",
    "object_lon",
    "object_date",
    "object_time"
  )

  if (!all(required_cols %in% colnames(eco_taxa_df))) {
    stop(
      paste(
        "The dataframe must contain the columns: ",
        paste(required_cols, collapse = ", ")
      )
    )
  }

  nrow_start <- nrow(eco_taxa_df)

  distinct_time_pos <- eco_taxa_df |>
    dplyr::distinct(
      object_lat,
      object_lon,
      object_date,
      object_time
    )

  # ensure the plan is reset when the function exits
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  # set up parallel backend
  future::plan(future::multisession, workers = workers)

  distinct_time_pos$is_day <- furrr::future_pmap_lgl(
    list(
      date = as.POSIXct(
        x = paste(
          distinct_time_pos$object_date,
          distinct_time_pos$object_time
        ),
        format = "%Y-%m-%d %H:%M:%S",
        tz = "UTC"
      ),
      lat = distinct_time_pos$object_lat,
      lon = distinct_time_pos$object_lon
    ),
    function(date, lat, lon) {
      SunCalcMeeus::is_daytime(
        date = date,
        geocode = tibble::tibble(
          lat = lat,
          lon = lon
        ),
        twilight = "nautical"
      )
    }
  )

  eco_taxa_df <- eco_taxa_df |>
    dplyr::left_join(
      y = distinct_time_pos,
      by = c(
        "object_lat",
        "object_lon",
        "object_date",
        "object_time"
      )
    ) |>
    dplyr::mutate(is_day = as.logical(is_day))

  # ensure that the number of rows in the input and output match
  agent <- pointblank::create_agent(tbl = eco_taxa_df) |>
    pointblank::col_vals_equal(
      columns = vars(nrow),
      value = nrow_start,
      preconditions = function(x) {
        x |>
          dplyr::summarize(nrow = dplyr::n())
      }
    ) |>
    pointblank::interrogate()

  if (!agent$validation_set$all_passed) {
    stop(
      "The number of rows in the input and output do not match."
    )
  }

  return(eco_taxa_df)

}