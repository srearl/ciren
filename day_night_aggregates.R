# day_night_aggregates <- eco_env |>
eco_env |>
  dplyr::group_by(
    object_date,
    object_time,
    d_n,
    is_day,
    object_lat,
    object_lon
  ) |>
  # dplyr::summarise(
  #   max(object_date),
  #   max(object_time),
  #   max(d_n),
  #   max(is_day),
  #   max(object_lat),
  #   max(object_lon)
  # ) |>
  # dplyr::ungroup() |>
  dplyr::mutate(
    sunrise = dplyr::if_else(
      !is.na(object_date) & !is.na(object_lat) & !is.na(object_lon),
      SunCalcMeeus::sunrise_time(
        date = as.Date(object_date),
        geocode = data.frame(
          lon = object_lon,
          lat = object_lat
        ),
        twilight = "nautical"
      ),
      NA
    ),
    sunset = dplyr::if_else(
      !is.na(object_date) & !is.na(object_lat) & !is.na(object_lon),
      SunCalcMeeus::sunset_time(
        date = as.Date(object_date),
        geocode = data.frame(
          lon = object_lon,
          lat = object_lat
        ),
        twilight = "nautical"
      ),
      NA
    )
  ) |>
  # dplyr::mutate(
  #   sunrise = SunCalcMeeus::sunrise_time(
  #     date = as.Date(object_date),
  #     geocode = data.frame(
  #       lon = object_lon,
  #       lat = object_lat
  #     ),
  #     twilight = "nautical"
  #   ),
  #   sunset = SunCalcMeeus::sunset_time(
  #     date = as.Date(object_date),
  #     geocode = data.frame(
  #       lon = object_lon,
  #       lat = object_lat
  #     ),
  #     twilight = "nautical"
  #   )
  # ) |>
  dplyr::select(-dplyr::matches("max"))


SunCalcMeeus::sunset_time(
  date     = as.Date("2022-04-01"),
  geocode  = data.frame(lon = -64, lat = 31.6),
  twilight = "nautical"
)

eco_env |>
dplyr::filter(
  grepl("nh", cruise, ignore.case = TRUE),
  grepl("m25", moc, ignore.case = TRUE)
  ) |> readr::write_csv("/tmp/beta.csv")


day_night_gradients <- eco_env |>
  dplyr::group_by(
    cruise,
    moc
  ) |>
  dplyr::summarise(
    max_object_date = max(object_date),
    max_object_time = max(object_time),
    max_d_n         = max(d_n),
    max_is_day      = max(is_day),
    max_object_lat  = max(object_lat, na.rm = TRUE),
    max_object_lon  = max(object_lon, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::rowwise() |>
  dplyr::mutate(
    sunrise = SunCalcMeeus::sunrise_time(
      date = as.Date(max_object_date),
      geocode = data.frame(
        lon = max_object_lon,
        lat = max_object_lat
      ),
      twilight = "nautical"
    ),
    sunset = SunCalcMeeus::sunset_time(
      date = as.Date(max_object_date),
      geocode = data.frame(
        lon = max_object_lon,
        lat = max_object_lat
      ),
      twilight = "nautical"
    )
  )
# dplyr::select(-dplyr::matches("max"))

SunCalcMeeus::sunset_time(
  date     = as.Date("2022-04-01"),
  geocode  = data.frame(lon = -64, lat = 31.6),
  twilight = "nautical"
)

alpha <- eco_env |>
  dplyr::group_by(
    cruise,
    moc
  ) |>
  dplyr::summarise(
    max_object_date = max(object_date),
    max_object_time = max(object_time),
    max_d_n         = max(d_n),
    max_is_day      = max(is_day),
    max_object_lat  = max(object_lat, na.rm = TRUE),
    max_object_lon  = max(object_lon, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::rowwise() |>
  dplyr::mutate(
    sunrise = SunCalcMeeus::sunrise_time(
      date = as.Date(max_object_date),
      geocode = data.frame(
        lon = max_object_lon,
        lat = max_object_lat
      ),
      twilight = "nautical"
    ),
    sunset = SunCalcMeeus::sunset_time(
      date = as.Date(max_object_date),
      geocode = data.frame(
        lon = max_object_lon,
        lat = max_object_lat
      ),
      twilight = "nautical"
    )
  )
# dplyr::select(-dplyr::matches("max"))
