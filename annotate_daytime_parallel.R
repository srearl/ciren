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
annotate_daytime_parallel <- function(
  eco_taxa_df,
  workers = future::availableCores() - 1
) {

  if (!all(c("cruise", "moc") %in% colnames(eco_taxa_df))) {
    stop("The dataframe must contain 'cruise' and 'moc' columns.")
  }

  nrow_start <- nrow(eco_taxa_df)

  # with join: 38.696 sec elapsed
  # sans join: 1676.436 sec elapsed

  grouped <- eco_taxa_df |>
    dplyr::group_by(
      cruise,
      moc
    ) |>
    dplyr::summarize(
      object_lat = mean(object_lat, na.rm = TRUE),
      object_lon = mean(object_lon, na.rm = TRUE),
      object_date = unique(object_date),
      object_time = unique(object_time),
      .groups = "drop"
    )

  # ensure that each group has only one unique date and time
  if (
    any(sapply(grouped$object_date, length) > 1) ||
      any(sapply(grouped$object_time, length) > 1)
  ) {
    stop("Each group must have only one unique date and time.")
  }

  # set up parallel backend
  future::plan(future::multisession, workers = workers)

  # apply `is_daytime` in parallel for each row
  grouped$is_day <- furrr::future_pmap_lgl(
    list(
      date = as.POSIXct(
        x = paste(grouped$object_date, grouped$object_time),
        format = "%Y-%m-%d %H:%M:%S",
        tz = "UTC"
      ),
      lat = grouped$object_lat,
      lon = grouped$object_lon
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
      y = grouped |>
        dplyr::select(
          cruise,
          moc,
          is_day
        ),
      by = c("cruise", "moc")
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(is_day = as.logical(is_day))

  # ensure that the number of rows in the input and output match
  agent <- pointblank::create_agent(tbl = eco_taxa_df) |>
    pointblank::col_vals_equal(
      columns = vars(nrow),
      value   = nrow_start,
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