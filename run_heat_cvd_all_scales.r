# ============================================
# Heat-associated CVD in Chicago: Mortality & ED
# End-to-end script (PHI-safe; no data embedded)
# Adds: county-level NO2 / PM2.5; MODIS NDVI at CA level
# Author: Peter Graffy
# Date: 2025-11-05
# ============================================

# ---- Packages ----
req <- c(
  "dplyr","readr","tidyr","stringr","purrr","ggplot2","lubridate",
  "mgcv","zoo","sf","tmap","here","glue","terra","boot","caret","gratia", "rlang", "forcats", "tibble"
)
to_install <- req[!req %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, Ncpus = max(1, parallel::detectCores()-1))
invisible(lapply(req, library, character.only=TRUE))

# ---- Paths & switches ----
root <- setwd("/Users/saborpete/Desktop/Peter/PhD/Dissertation/Analysis/Code and Data/Heat_CVD_Chicago")


# --- Robust file finder (searches in project root) ---
find_file <- function(candidates,
                      root_dir = here::here(),
                      extra_dirs = NULL,
                      verbose = TRUE) {
  search_dirs <- unique(c(root_dir, getwd(), extra_dirs))
  # 1) exact filename checks
  exact_paths <- unlist(lapply(search_dirs, function(d) file.path(d, candidates)))
  hit_exact <- exact_paths[file.exists(exact_paths)]
  if (length(hit_exact)) {
    p <- normalizePath(hit_exact[1], mustWork = TRUE)
    if (verbose) message("Found file: ", p)
    return(p)
  }
  # 2) fallback: case-insensitive single-level pattern search
  for (cand in candidates) {
    for (d in search_dirs) {
      files <- list.files(d, pattern = paste0("^", gsub("\\.", "\\\\.", cand), "$"),
                          ignore.case = TRUE, full.names = TRUE, recursive = FALSE)
      if (length(files)) {
        p <- normalizePath(files[1], mustWork = TRUE)
        if (verbose) message("Found (ci) file: ", p)
        return(p)
      }
    }
  }
  if (verbose) {
    message("Could not find any of:\n  ",
            paste(unique(file.path(search_dirs, candidates)), collapse = "\n  "))
  }
  return(NULL)
}

data_dir <- root

out_fig <- file.path(root, "output","figures"); dir.create(out_fig, recursive=TRUE, showWarnings=FALSE)
out_tab <- file.path(root, "output","tables");  dir.create(out_tab, recursive=TRUE, showWarnings=FALSE)
out_log <- file.path(root, "output","logs");    dir.create(out_log, recursive=TRUE, showWarnings=FALSE)
out_tmp <- file.path(root, "output","scratch"); dir.create(out_tmp, recursive=TRUE, showWarnings=FALSE)

paths <- list(
  deaths   = find_file(c("deaths_outcomes_complete.csv")),
  ed       = find_file(c("ed_outcomes_complete.csv")),
  acs      = find_file(c("chicago_acs_2010-2022.csv","acs_ca_2010-2022.csv")),
  xwalk    = find_file(c("areas_tracts2010_joined.csv","areas_tracts2010_joined.csv")),
  temps    = find_file(c("CA_temps_rh_90-23.csv")),
  cas      = find_file(c("comm_areas.geojson")),
  # Air pollution
  no2_old  = find_file(c("conus_county_no2_2005_2020.csv")),
  no2_new  = find_file(c("no2_county_year.csv")),
  pm25     = find_file(c("conus_county_pm25_1998_2023.csv","pm25_county_year.csv")),
  # NDVI
  ndvi_year_csv  = find_file(c("ndvi_year_ca.csv")),
  ndvi_month_csv = find_file(c("ndvi_month_ca.csv"))
)

cfg <- list(
  start_year = 2011, end_year = 2022,
  warm_months_only = TRUE,               # May–September
  Covid_indicator = TRUE,
  Covid_exclusion = TRUE,
  adj_air_green   = TRUE,                # include NO2/PM2.5/NDVI if built
  # Spline controls
  k_tmax_year  = 5, k_tmax_month = 6, k_tmax_day = 7,
  k_age = 5, k_hum = 6, k_time_year = 5, k_time_date = 10,
  # Threshold sequences
  thr_year  = seq(18, 35, by = 0.05),
  thr_month = seq(18, 35, by = 0.25),
  thr_day   = seq(10, 42, by = 0.5),
  # Bootstrap
  # Bootstrap controls (full grid; no down-sampling)
  boot_R = 200,
  ci_workers = max(1, parallel::detectCores() - 1),
  boot_parallel = if (.Platform$OS.type == "windows") "snow" else "multicore",
  ci_seed = 20251105,
  # NDVI options (only used if downloading)
  ndvi_download = FALSE,
  ndvi_product  = "MOD13A2",            # 1km 16-day NDVI
  ndvi_start    = "2011-05-01",
  ndvi_end      = "2022-09-30",
  ndvi_user     = Sys.getenv("EARTHDATA_USER"),
  ndvi_pass     = Sys.getenv("EARTHDATA_PASS")
)

expected_flags <- c("CVD","CHD","MI","Stroke","HBP","Diabetes")
COOK_GEOID <- "17031"

# ---- Helpers ----
try_read_csv <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  out <- try(suppressWarnings(suppressMessages(
    readr::read_csv(path, show_col_types = FALSE)
  )), silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

diverging_breaks <- function(x, n_breaks = 11) {
  neg <- x[x < 0]; pos <- x[x > 0]
  n_side <- floor((n_breaks - 1)/2)
  nb <- if (length(neg)) quantile(neg, probs=seq(0,1,length.out=n_side+1), na.rm=TRUE) else 0
  pb <- if (length(pos)) quantile(pos, probs=seq(0,1,length.out=n_side+1), na.rm=TRUE) else 0
  unique(c(nb[-length(nb)], 0, pb[-1]))
}

safe_diverging_scale <- function(x, n_ticks = 9, min_span = 1e-6) {
  xf <- x[is.finite(x)]
  if (!length(xf)) {
    return(list(limits = c(-1, 1), ticks = c(-0.5, 0, 0.5), ok = FALSE))
  }
  
  rng <- range(xf, na.rm = TRUE)
  a <- max(abs(rng))
  if (!is.finite(a) || a < min_span) a <- 1
  
  lims <- c(-a, a)
  
  # candidate ticks
  ticks <- pretty(lims, n = n_ticks)
  
  # keep ticks strictly inside limits (tmap complains otherwise)
  eps <- .Machine$double.eps^0.25
  ticks <- ticks[ticks > (lims[1] + eps) & ticks < (lims[2] - eps)]
  
  # ensure 0 is shown if it lies inside the limits
  if (lims[1] < 0 && 0 < lims[2]) ticks <- sort(unique(c(ticks, 0)))
  
  # if clamping left too few ticks, build a simple interior sequence
  if (length(ticks) < 3) {
    ticks <- seq(lims[1] + (a/4), lims[2] - (a/4), length.out = 3)
    # include 0 when possible
    if (lims[1] < 0 && 0 < lims[2]) ticks[2] <- 0
  }
  
  list(limits = lims, ticks = ticks, ok = TRUE)
}

mk_rates <- function(df) {
  df %>%
    mutate(
      CVD_per_100k      = (CVD_count / total_pop) * 1e5,
      CHD_per_100k      = (CHD_count / total_pop) * 1e5,
      MI_per_100k       = (MI_count  / total_pop) * 1e5,
      Stroke_per_100k   = (Stroke_count / total_pop) * 1e5,
      HBP_per_100k      = (HBP_count / total_pop) * 1e5,
      Diabetes_per_100k = (Diabetes_count / total_pop) * 1e5,
      Deaths_per_100k   = (all_deaths / total_pop) * 1e5,
      CVD_deaths_ratio  = ifelse(all_deaths>0, CVD_count / all_deaths, NA_real_)
    )
}

# ---- Safe scaler for added covariates (avoid NaN when sd=0 or all NA) ----
.safe_scale <- function(x) {
  if (all(is.na(x))) return(x)
  sdx <- stats::sd(x, na.rm = TRUE)
  if (is.na(sdx) || sdx == 0) return(x - mean(x, na.rm = TRUE))
  as.numeric((x - mean(x, na.rm = TRUE)) / sdx)
}

z_air_green <- function(df) {
  cols <- intersect(c("ndvi","no2","pm25"), names(df))
  if (!length(cols)) return(df)
  dplyr::mutate(df, dplyr::across(dplyr::all_of(cols), .safe_scale))
}

# ---------------------------
# Bootstrap CI for threshold scans using predictions or observed
# ---------------------------
boot_excess_stat <- function(data, indices, temp_var, pred_var, group = "community", threshold) {
  d <- data[indices, ]
  d <- d %>% mutate(heat = .data[[temp_var]] > threshold)
  d %>%
    group_by(.data[[group]]) %>%
    summarise(
      mean_hot  = mean(.data[[pred_var]][heat],  na.rm=TRUE),
      mean_cool = mean(.data[[pred_var]][!heat], na.rm=TRUE),
      diff = mean_hot - mean_cool,
      .groups="drop"
    ) %>%
    summarise(avg_excess = mean(diff, na.rm=TRUE), .groups="drop") %>%
    pull(avg_excess)
}

# ---- Helper: linear interpolation of CI columns from coarse -> full grid ----
interpolate_ci_to_full <- function(ci_coarse, full_thresholds) {
  if (is.null(ci_coarse) || nrow(ci_coarse) < 2) {
    return(tibble::tibble(threshold = full_thresholds,
                          lower_ci = NA_real_, upper_ci = NA_real_))
  }
  f <- function(y) stats::approx(ci_coarse$threshold, y,
                                 xout = full_thresholds,
                                 method = "linear", rule = 2)$y
  tibble::tibble(threshold = full_thresholds,
                 lower_ci = f(ci_coarse$lower_ci),
                 upper_ci = f(ci_coarse$upper_ci))
}

# ---- FAST bootstrap CIs: stride + parallel + interpolation ----
bootstrap_threshold_curve <- function(df, temp_var, pred_var, thr_seq, R = 200, group = "community") {
  n <- length(thr_seq)
  pb_local <- utils::txtProgressBar(min = 0, max = n, style = 3)
  on.exit(try(close(pb_local), silent = TRUE), add = TRUE)
  set.seed(cfg$ci_seed %||% 1L)
  
  res <- vector("list", n)
  for (k in seq_along(thr_seq)) {
    th <- thr_seq[k]
    b <- boot::boot(
      data = df,
      statistic = function(d, i) boot_excess_stat(
        d, i, temp_var = temp_var, pred_var = pred_var, group = group, threshold = th
      ),
      R = R,
      parallel = cfg$boot_parallel %||% "no",
      ncpus    = cfg$ci_workers %||% 1
    )
    ci <- try(boot::boot.ci(b, type = "perc")$percent[4:5], silent = TRUE)
    res[[k]] <- tibble::tibble(
      threshold = th,
      lower_ci  = if (inherits(ci,"try-error")) NA_real_ else ci[1],
      upper_ci  = if (inherits(ci,"try-error")) NA_real_ else ci[2]
    )
    utils::setTxtProgressBar(pb_local, k)
  }
  dplyr::bind_rows(res)
}

# ---------------------------
# Scan thresholds (observed and/or predicted)
# ---------------------------
scan_threshold <- function(df, temp_var, rate_var, thr_seq, pred_vec = NULL, group = "community") {
  y_in <- if (is.null(pred_vec)) df[[rate_var]] else pred_vec
  res <- purrr::map_df(thr_seq, function(th) {
    df %>%
      mutate(heat = .data[[temp_var]] > th, y = y_in) %>%
      group_by(.data[[group]]) %>%
      summarise(
        mean_hot  = mean(y[heat],  na.rm=TRUE),
        mean_cool = mean(y[!heat], na.rm=TRUE),
        diff = mean_hot - mean_cool, .groups="drop"
      ) %>%
      summarise(threshold = th,
                avg_excess_rate = mean(diff, na.rm=TRUE),
                sd_excess_rate  = sd(diff, na.rm=TRUE),
                n = dplyr::n())
  })
  peak <- res %>% filter(avg_excess_rate == max(avg_excess_rate, na.rm=TRUE))
  highest_pos <- res %>% filter(avg_excess_rate > 0) %>% arrange(desc(threshold)) %>% slice_head(n=1)
  list(grid=res, peak=peak, highest_positive=highest_pos)
}

# ---------------------------
# Flexible threshold plotter (adds ribbon if CI cols provided)
# ---------------------------
# ---- Threshold plot with tight CI + dotted max line + legend values at bottom ----
plot_thr <- function(grid_df,
                     peak,                 # tibble with $threshold
                     highest_pos = NULL,   # tibble with $threshold (optional)
                     xlab,
                     file,
                     ci_lo = NULL, ci_hi = NULL,  # optional CI column names
                     subtitle = "Chicago, 2011–2022") {
  
  stopifnot(all(c("threshold","avg_excess_rate") %in% names(grid_df)))
  
  # Optional ribbon only if both CI columns exist
  has_ci <- !is.null(ci_lo) && !is.null(ci_hi) &&
    all(c(ci_lo, ci_hi) %in% names(grid_df))
  
  # Format threshold values for legend labels (show 1 decimal; fall back to NA if missing)
  fmt <- function(x) if (length(x) && is.finite(x[1])) sprintf("%.1f°C", x[1]) else "NA"
  peak_val <- fmt(peak$threshold)
  high_val <- fmt(if (is.null(highest_pos)) NA_real_ else highest_pos$threshold)
  
  # Build legend labels *with values*
  peak_lab <- paste0("Peak threshold (max excess): ", peak_val)
  high_lab <- paste0("Highest threshold with > 0 excess: ", high_val)
  
  vlines <- tibble::tibble(
    x     = c(peak$threshold[1], if (!is.null(highest_pos)) highest_pos$threshold[1] else NA_real_),
    label = factor(c(peak_lab, high_lab), levels = c(peak_lab, high_lab))
  ) %>% dplyr::filter(is.finite(x))
  
  g <- ggplot(grid_df, aes(threshold, avg_excess_rate)) +
    # CI ribbon drawn *under* the line, subtle alpha so it hugs visually
    { if (has_ci) geom_ribbon(aes(ymin = .data[[ci_lo]], ymax = .data[[ci_hi]]),
                              alpha = 0.12) } +
    geom_line(linewidth = 1.0, na.rm = TRUE) +  # line on top so the ribbon appears tight
    # Two vertical lines with distinct linetypes (max = dotted, highest-pos = dashed)
    geom_vline(
      data = vlines,
      aes(xintercept = x, linetype = label),
      linewidth = 0.9,
      show.legend = TRUE
    ) +
    scale_linetype_manual(
      name = NULL,
      values = setNames(c("dotted", "dashed"), c(peak_lab, high_lab))  # dotted = “fine dots”
    ) +
    labs(x = xlab,
         y = "Average heat-associated excess per 100k",
         subtitle = subtitle) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",           # legend at bottom
      legend.box = "vertical",
      legend.margin = margin(t = 4, r = 4, b = 4, l = 4)
    )
  
  # keep x-range tight to scanned range using finite points only
  xr <- range(grid_df$threshold[is.finite(grid_df$avg_excess_rate)], na.rm = TRUE)
  g <- g + coord_cartesian(xlim = xr)
  
  ggsave(file, g, width = 7.5, height = 5.5, dpi = 300)
  invisible(g)
}


# ---- Air pollution loaders (county-level → annual; Cook only) ----
load_no2_annual <- function() {
  new <- if (!is.null(paths$no2_new)) readr::read_csv(paths$no2_new, show_col_types = FALSE) else NULL
  old <- if (!is.null(paths$no2_old)) readr::read_csv(paths$no2_old, show_col_types = FALSE) else NULL
  if (is.null(new) && is.null(old)) return(NULL)
  norm <- function(df) {
    nm <- names(df)
    if (all(c("GEOID","year","no2_mean") %in% nm)) {
      return(df %>% transmute(year = as.integer(year),
                              GEOID = sprintf("%05s", gsub("\\.0+$","", as.character(GEOID))),
                              no2 = as.numeric(no2_mean)))
    }
    if (all(c("year","GEOID","mean_no2") %in% nm)) {
      return(df %>% transmute(year = as.integer(year),
                              GEOID = sprintf("%05s", gsub("\\.0+$","", as.character(GEOID))),
                              no2 = as.numeric(mean_no2)))
    }
    stop("NO2: Unrecognized column layout: ", paste(nm, collapse=", "))
  }
  all <- dplyr::bind_rows(
    if (!is.null(old)) tryCatch(norm(old), error = function(e) NULL),
    if (!is.null(new)) tryCatch(norm(new), error = function(e) NULL)
  )
  if (is.null(all) || !nrow(all)) return(NULL)
  all %>%
    filter(GEOID == COOK_GEOID) %>%
    arrange(year) %>%
    group_by(year) %>% summarise(no2 = dplyr::last(no2), .groups="drop") %>%
    filter(year >= cfg$start_year, year <= cfg$end_year)
}

load_pm25_annual <- function() {
  df <- if (!is.null(paths$pm25)) readr::read_csv(paths$pm25, show_col_types = FALSE) else NULL
  if (is.null(df)) return(NULL)
  nm <- names(df)
  val_col <- if ("mean_pm25" %in% nm) "mean_pm25" else if ("pm25_mean" %in% nm) "pm25_mean" else NULL
  if (is.null(val_col)) stop("PM2.5: could not find mean_pm25/pm25_mean in: ", paste(nm, collapse=", "))
  df %>%
    transmute(year = as.integer(year),
              GEOID = sprintf("%05s", gsub("\\.0+$","", as.character(GEOID))),
              pm25  = as.numeric(.data[[val_col]])) %>%
    filter(GEOID == COOK_GEOID,
           year >= cfg$start_year, year <= cfg$end_year) %>%
    arrange(year)
}

# ---- NDVI builders (use CSVs if present; otherwise skip) ----
have_ndvi_year  <- function() !is.null(paths$ndvi_year_csv)  && file.exists(paths$ndvi_year_csv)
have_ndvi_month <- function() !is.null(paths$ndvi_month_csv) && file.exists(paths$ndvi_month_csv)

normalize_ndvi_year <- function(df) {
  # Accept common shapes and return community, year, ndvi (mean)
  nm <- names(df)
  if (all(c("community","year","ndvi_mean") %in% nm)) {
    return(df %>% transmute(community, year = as.integer(year), ndvi = as.numeric(ndvi_mean)))
  }
  if (all(c("community","year","ndvi") %in% nm)) {
    return(df %>% transmute(community, year = as.integer(year), ndvi = as.numeric(ndvi)))
  }
  stop("NDVI (year) CSV lacks expected columns (community/year/ndvi_mean or ndvi).")
}

build_ndvi_tables <- function() {
  if (have_ndvi_year() || have_ndvi_month()) {
    message("Using existing NDVI CSV(s).")
    ndvi_year_raw  <- if (have_ndvi_year())  readr::read_csv(paths$ndvi_year_csv,  show_col_types = FALSE) else NULL
    ndvi_month_raw <- if (have_ndvi_month()) readr::read_csv(paths$ndvi_month_csv, show_col_types = FALSE) else NULL
    ndvi_year  <- if (!is.null(ndvi_year_raw))  normalize_ndvi_year(ndvi_year_raw) else NULL
    return(list(year = ndvi_year, month = ndvi_month_raw))
  }
  if (!cfg$ndvi_download) {
    message("NDVI download disabled and no CSVs found — proceeding without NDVI.")
    return(NULL)
  }
  stop("NDVI download requested but EARTHDATA credentials are not set.")
}

# =================
# Mapping (tmap v4)
# =================
# =================
# Mapping (tmap v4): red for +, blue for -
# =================
make_rate_map <- function(cas_sf,
                          panel_df,
                          threshold,
                          legend_title,
                          out_file,
                          rate_col = "pred",
                          temp_col = "mean_tmax") {
  
  stopifnot(all(c("community", temp_col, rate_col) %in% names(panel_df)))
  
  map_dat <- panel_df %>%
    dplyr::group_by(community) %>%
    dplyr::mutate(is_hot = .data[[temp_col]] > threshold) %>%
    dplyr::summarise(
      mean_hot  = mean(.data[[rate_col]][is_hot],  na.rm = TRUE),
      mean_cool = mean(.data[[rate_col]][!is_hot], na.rm = TRUE),
      excess    = mean_hot - mean_cool,
      .groups = "drop"
    ) %>%
    dplyr::right_join(cas_sf, by = "community") %>%
    sf::st_as_sf()
  
  # Diverging breaks & limits (safe + symmetric around 0)
  scale_info <- safe_diverging_scale(map_dat$excess)
  
  # If nothing finite, skip gracefully
  if (!scale_info$ok) {
    message("make_rate_map: no finite 'excess' values — map skipped.")
    return(invisible(NULL))
  }
  
  # RdBu reversed: negatives → blue, positives → red
  p <- tm_shape(map_dat) +
    tm_polygons(
      fill = "excess",
      fill.scale = tm_scale_continuous(
        values   = "-brewer.RdBu",   # blue for −, red for +
        midpoint = 0,
        limits   = scale_info$limits,
        ticks    = scale_info$ticks
      ),
      fill.legend = tm_legend(
        title = legend_title,
        reverse = FALSE
      ),
      col = "white",
      lwd = 0.4,
      col_alpha = 0.8
    ) +
    tm_layout(
      frame = FALSE,
      legend.outside = TRUE
    )
  
  tmap::tmap_save(p, filename = out_file, width = 1400, height = 1000, dpi = 200)
}


# ---- Core resources ----
acs    <- readr::read_csv(paths$acs,   show_col_types = FALSE)
xwalk  <- readr::read_csv(paths$xwalk, show_col_types = FALSE)
temps  <- readr::read_csv(paths$temps, show_col_types = FALSE)
cas    <- sf::st_read(paths$cas, quiet = TRUE)
stopifnot(!is.null(acs), !is.null(xwalk), !is.null(temps), !inherits(cas,"try-error"))

# ----- ACS summarize (community × year) -----
acs_sel <- acs[!(acs$variable %in% c("Total Pop","Hispanic origin")), ]

acs_join <- dplyr::left_join(acs_sel, xwalk, by = c("GEOID" = "geoid10"))
acs_summary <- acs_join %>%
  group_by(community, year) %>%
  summarise(
    total_pop     = sum(estimate[variable == "Total population"], na.rm = TRUE),
    mean_college  = mean(percent[variable == "Bachelor's degree"], na.rm = TRUE),
    mean_hs       = mean(percent[variable == "High school"], na.rm = TRUE),
    mean_assoc    = mean(percent[variable == "Associate's degree"], na.rm = TRUE),
    mean_white    = mean(percent[variable == "White"], na.rm = TRUE),
    mean_black    = mean(percent[variable == "Black"], na.rm = TRUE),
    mean_asian    = mean(percent[variable == "Asian"], na.rm = TRUE),
    mean_hispanic = mean(percent[variable == "Hispanic"], na.rm = TRUE),
    mean_employed = mean(percent[variable == "Employed"], na.rm = TRUE),
    mean_unemployed = mean(percent[variable == "Unemployed"], na.rm = TRUE),
    mean_male     = mean(percent[variable == "Male"], na.rm = TRUE),
    mean_female   = mean(percent[variable == "Female"], na.rm = TRUE),
    median_age    = median(estimate[variable == "Median age"], na.rm = TRUE),
    median_income = median(estimate[variable == "Median income"], na.rm = TRUE),
    .groups = "drop"
  )

# Temps → panels
temps <- temps %>% mutate(date = as.Date(date))
if (cfg$warm_months_only) temps <- temps %>% filter(lubridate::month(date) %in% 5:9)
temps <- temps %>% filter(year(date) >= cfg$start_year, year(date) <= cfg$end_year)
temps_daily <- temps %>% mutate(year=year(date), month=month(date), day=day(date))
temps_monthly <- temps_daily %>%
  group_by(community, year, month) %>%
  summarise(across(c(tmax,tmin,tmean,vp,humidity), ~mean(.x, na.rm=TRUE), .names="mean_{.col}"), .groups="drop")
temps_yearly <- temps_daily %>%
  group_by(community, year) %>%
  summarise(across(c(tmax,tmin,tmean,vp,humidity), ~mean(.x, na.rm=TRUE), .names="mean_{.col}"), .groups="drop")

# Air pollution tables (Cook County → join by year)
no2_tbl  <- load_no2_annual()  # year, no2
pm25_tbl <- load_pm25_annual() # year, pm25
message("NO2 rows: ", if (is.null(no2_tbl)) 0 else nrow(no2_tbl),
        " | PM2.5 rows: ", if (is.null(pm25_tbl)) 0 else nrow(pm25_tbl))

# NDVI tables (community-year)
ndvi_tables <- if (cfg$adj_air_green) build_ndvi_tables() else NULL
ndvi_year   <- if (!is.null(ndvi_tables)) ndvi_tables$year  else NULL

# ======================================================
# Pipeline for a given outcome: "mortality" | "ed"
# ======================================================
# ---- Progress bar (txt) ----
pb <- local({
  bar <- NULL; total <- 0; i <- 0
  list(
    start = function(n, header = NULL) {
      total <<- n; i <<- 0
      if (!is.null(header)) message(header)
      bar <<- utils::txtProgressBar(min = 0, max = n, style = 3)
    },
    tick = function(msg = NULL) {
      i <<- i + 1
      if (!is.null(msg)) message(sprintf("[%03d/%03d] %s", i, total, msg))
      if (!is.null(bar)) utils::setTxtProgressBar(bar, i)
    },
    end = function() {
      if (!is.null(bar)) close(bar)
      bar <<- NULL
    }
  )
})

try_read_csv <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  out <- try(suppressMessages(readr::read_csv(path, show_col_types = FALSE)), silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

safe_rename_response <- function(df, rate_col) {
  # If we're modeling CVD itself, do nothing.
  if (identical(rate_col, "CVD_per_100k")) return(df)
  
  if (!rate_col %in% names(df)) {
    stop("rate_col not found in data: ", rate_col,
         " | available: ", paste(names(df), collapse = ", "))
  }
  
  # Remove any existing CVD_per_100k to avoid duplicates, then rename selected rate to that name
  df %>%
    dplyr::select(-dplyr::any_of("CVD_per_100k")) %>%
    dplyr::rename(CVD_per_100k = !! rlang::sym(rate_col))
}

# --- Make/derive a 'community' column from many possibilities ---
ensure_community <- function(df) {
  nm <- names(df)
  
  # already present?
  if ("community" %in% nm) return(df)
  
  # common name columns
  name_cands <- intersect(nm, c("CA_NAME","ca_name","CAName",
                                "community_area_name","community_area","Community","COMMUNITY"))
  if (length(name_cands)) {
    return(dplyr::rename(df, community = !! rlang::sym(name_cands[1])))
  }
  
  # community area numeric → join to cas (expects cas has area_num_1 and community)
  num_cands <- intersect(nm, c("COMMUNITY_AREA","community_area_number","area_num_1",
                               "area_numbe","ca_num","CA_NUM"))
  if (length(num_cands) && "area_num_1" %in% names(cas)) {
    lut <- sf::st_drop_geometry(cas) %>%
      dplyr::select(area_num_1, community) %>%
      dplyr::distinct()
    return(df %>%
             dplyr::left_join(lut, by = setNames("area_num_1", num_cands[1])))
  }
  
  # tract GEOID → join to xwalk (expects xwalk has geoid10 and community)
  geoid_cands <- intersect(nm, c("GEOID","geoID","geoid10","TRACT_GEOID","tract_geoid"))
  if (length(geoid_cands) && "geoid10" %in% names(xwalk)) {
    gx <- xwalk %>% dplyr::select(geoid10, community) %>% dplyr::distinct()
    # normalize to 11-digit tract id if it looks longer (e.g., block group)
    col <- geoid_cands[1]
    df2 <- df %>% dplyr::mutate(
      ..geoid_tmp.. = gsub("[^0-9]", "", as.character(.data[[col]])),
      ..geoid_tmp.. = substr(..geoid_tmp.., 1, 11)
    )
    out <- df2 %>% dplyr::left_join(gx, by = c("..geoid_tmp.." = "geoid10")) %>%
      dplyr::select(-..geoid_tmp..)
    return(out)
  }
  
  stop("Could not derive 'community'. Provide CA_NAME, a community-area code, or tract GEOID.")
}

# --- Where to save things (adjust if you already have out_* elsewhere) ---
out_dir  <- file.path("output")                  # or cfg$out_dir
out_mod  <- file.path(out_dir, "models")
out_tab  <- file.path(out_dir, "tables")
out_fig  <- file.path(out_dir, "figures")
out_log  <- file.path(out_dir, "logs")
dir.create(out_mod, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log, recursive = TRUE, showWarnings = FALSE)

# --- Save a GAM with covariate baselines for later partial-effect plots ---
save_gam_with_means <- function(fit, path) {
  mf <- fit$model
  num_means <- vapply(mf, function(x) if (is.numeric(x)) mean(x, na.rm=TRUE) else NA_real_, 0)
  fac_bases <- lapply(mf, function(x) if (is.factor(x)) levels(x)[1] else NULL)
  obj <- list(
    fit = fit,
    num_means = num_means[!is.na(num_means)],
    fac_bases = fac_bases,
    linkinv   = fit$family$linkinv
  )
  saveRDS(obj, path)
  invisible(path)
}

# --- Extract edf / p for the temperature smooth + adj.R² + dev.expl ---
summarize_gam <- function(fit, temp_pattern) {
  sm <- summary(fit)
  st <- sm$s.table
  row <- st[grep(temp_pattern, rownames(st))[1], , drop = FALSE]
  data.frame(
    edf   = if (nrow(row)) unname(row[,"edf"]) else NA_real_,
    p     = if (nrow(row)) unname(row[,"p-value"]) else NA_real_,
    r2adj = sm$r.sq,
    dev   = sm$dev.expl
  )
}

# --- Append or create a CSV safely (no duplicated headers) ---
append_csv_row <- function(df_row, path) {
  if (!file.exists(path)) {
    readr::write_csv(df_row, path)
  } else {
    readr::write_csv(df_row, path, append = TRUE)
  }
}



run_pipeline <- function(outcome_type = c("mortality","ed")) {
  outcome_type <- match.arg(outcome_type)
  message(paste0("▶ Running pipeline: ", toupper(outcome_type)))
  
  # load + date
  if (outcome_type=="mortality") {
    df <- try_read_csv(paths$deaths); stopifnot(!is.null(df))
    date_col <- c("DOD","DATE")[c("DOD","DATE") %in% names(df)][1]
  } else {
    df <- try_read_csv(paths$ed); stopifnot(!is.null(df))
    date_col <- c("ENC_ADMIT_DATE","ED_DATE","DATE")[c("ENC_ADMIT_DATE","ED_DATE","DATE") %in% names(df)][1]
    if ("geoID" %in% names(df) && !"GEOID" %in% names(df)) df <- dplyr::rename(df, GEOID = geoID)
  }
  stopifnot(!is.null(date_col))
  df <- df %>% dplyr::mutate(DATE = as.Date(.data[[date_col]]))
  
  # ensure 'community'
  df <- ensure_community(df)
  if (!"community" %in% names(df)) stop("Outcome file: could not make 'community' column.")
  
  # flags present?
  present <- intersect(expected_flags, names(df))
  if (!length(present)) stop("No CVD/CHD/MI/Stroke flags found in the outcome file.")
  
  # time filters
  if (cfg$warm_months_only) df <- df %>% dplyr::filter(lubridate::month(DATE) %in% 5:9)
  df <- df %>% dplyr::filter(lubridate::year(DATE) >= cfg$start_year,
                             lubridate::year(DATE) <= cfg$end_year)
    # ---- Aggregate daily / monthly / yearly ----
  message("• Building daily/monthly/yearly panels …")
  # ---- Aggregate daily / monthly / yearly ----
  counts_daily <- df %>%
    dplyr::group_by(community, DATE) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(present), ~sum(.x, na.rm=TRUE), .names="{.col}_count"),
      all_deaths = dplyr::n(), .groups="drop"
    ) %>%
    dplyr::mutate(year = lubridate::year(DATE),
                  month = lubridate::month(DATE),
                  day = lubridate::day(DATE))
  
  counts_monthly <- counts_daily %>%
    dplyr::group_by(community, year, month) %>%
    dplyr::summarise(dplyr::across(ends_with("_count"), sum, na.rm=TRUE),
                     all_deaths = sum(all_deaths, na.rm=TRUE), .groups="drop")
  
  counts_yearly <- counts_daily %>%
    dplyr::group_by(community, year) %>%
    dplyr::summarise(dplyr::across(ends_with("_count"), sum, na.rm=TRUE),
                     all_deaths = sum(all_deaths, na.rm=TRUE), .groups="drop")
  
  panel_yearly <- acs_summary %>%
    dplyr::inner_join(counts_yearly, by=c("community","year")) %>%
    dplyr::inner_join(temps_yearly,  by=c("community","year")) %>% mk_rates()
  
  panel_monthly <- acs_summary %>%
    dplyr::inner_join(counts_monthly, by=c("community","year")) %>%
    dplyr::inner_join(temps_monthly,  by=c("community","year","month")) %>%
    dplyr::mutate(date_num = as.numeric(lubridate::ymd(paste(year, month, "01")))) %>% mk_rates()
  
  panel_daily <- temps_daily %>%
    dplyr::left_join(acs_summary,  by=c("community","year")) %>%
    dplyr::left_join(counts_daily, by=c("community","year","month","day")) %>%
    dplyr::select(-DATE) %>%
    dplyr::mutate(dplyr::across(dplyr::ends_with("_count"),
                                ~ tidyr::replace_na(.x, 0))) %>%
    mk_rates() %>%
    dplyr::mutate(date_num = as.numeric(date)) %>%
    dplyr::arrange(community, date) %>%
    dplyr::group_by(community) %>%
    dplyr::mutate(tmax_3day_avg = zoo::rollapply(tmax, width=3, FUN=mean, align="right", fill=NA, na.rm=TRUE)) %>%
    dplyr::ungroup()
  
  # ---- Add NO2/PM2.5/NDVI + standardize ----
  if (!is.null(pm25_tbl)) { panel_yearly <- panel_yearly %>% dplyr::left_join(pm25_tbl, by="year")
  panel_monthly<- panel_monthly%>% dplyr::left_join(pm25_tbl, by="year") }
  if (!is.null(no2_tbl))  { panel_yearly <- panel_yearly %>% dplyr::left_join(no2_tbl,  by="year")
  panel_monthly<- panel_monthly%>% dplyr::left_join(no2_tbl,  by="year") }
  if (!is.null(ndvi_year)) {
    panel_yearly  <- panel_yearly  %>% dplyr::left_join(ndvi_year,  by=c("community","year"))
    panel_monthly <- panel_monthly %>% dplyr::left_join(ndvi_year,  by=c("community","year"))
    panel_daily   <- panel_daily   %>% dplyr::left_join(ndvi_year,  by=c("community","year"))
  }
  panel_yearly  <- z_air_green(panel_yearly)
  panel_monthly <- z_air_green(panel_monthly)
  panel_daily   <- z_air_green(panel_daily)
  
  if (cfg$Covid_indicator) {
    panel_yearly  <- panel_yearly  %>% dplyr::mutate(covid_period = year >= 2020)
    panel_monthly <- panel_monthly %>% dplyr::mutate(covid_period = year >= 2020)
    panel_daily   <- panel_daily   %>% dplyr::mutate(covid_period = year >= 2020)
  }
  
  # ---- Endpoint loop with progress ----
  rate_cols <- c("CVD_per_100k","CHD_per_100k","MI_per_100k","Stroke_per_100k")
  prefix <- ifelse(outcome_type=="mortality","MORT","ED")
  
  # ~31 ticks per endpoint (fit 4, scans 8, CIs 8, plots 8, map 1, sens 2-ish)
  TICKS_PER_EP <- 31
  pb$total <- NULL
  pb$start(length(rate_cols) * TICKS_PER_EP,
           header = paste0("• Processing endpoints for ", toupper(outcome_type), " …"))
  
  for (rate_col in rate_cols) {
    if (!rate_col %in% names(panel_yearly)) { next }
    
    short <- sub("_per_100k$","", rate_col)
    message(sprintf("  → %s (%s)", short, toupper(outcome_type)))
    
    py  <- safe_rename_response(panel_yearly,  rate_col)
    pm  <- safe_rename_response(panel_monthly, rate_col)
    pd  <- safe_rename_response(panel_daily,   rate_col)
    pdma <- pd %>% dplyr::filter(!is.na(tmax_3day_avg))
    
    # Model defs
    gam_yearly <- function(df) {
      extra <- intersect(c("pm25","no2","ndvi"), names(df))
      rhs <- paste(
        sprintf("s(mean_tmax, k=%d)", cfg$k_tmax_year),
        sprintf("s(median_age, k=%d)", cfg$k_age),
        sprintf("s(mean_humidity, k=%d)", cfg$k_hum),
        "median_income", "mean_college", "mean_hs",
        "mean_white", "mean_black", "mean_hispanic", "mean_asian",
        "mean_employed", "mean_unemployed", "mean_male", "mean_female",
        sprintf("s(year, bs='cs', k=%d)", cfg$k_time_year),
        if (length(extra)) paste(extra, collapse = " + ") else NULL,
        sep = " + "
      )
      mgcv::gam(as.formula(paste0("CVD_per_100k ~ ", rhs)),
                data=df, family = mgcv::nb(link="log"), method="REML", select=TRUE)
    }
    gam_monthly <- function(df) {
      extra <- intersect(c("pm25","no2","ndvi"), names(df))
      rhs <- paste(
        sprintf("s(mean_tmax, k=%d)", cfg$k_tmax_month),
        sprintf("s(median_age, k=%d)", cfg$k_age),
        sprintf("s(mean_humidity, k=%d)", cfg$k_hum),
        "median_income", "mean_college", "mean_hs",
        "mean_white", "mean_black", "mean_hispanic", "mean_asian",
        "mean_employed", "mean_unemployed", "mean_male", "mean_female",
        sprintf("s(date_num, bs='cs', k=%d)", cfg$k_time_date),
        if (length(extra)) paste(extra, collapse = " + ") else NULL,
        sep = " + "
      )
      mgcv::gam(as.formula(paste0("CVD_per_100k ~ ", rhs)),
                data=df, family = mgcv::nb(link="log"), method="REML", select=TRUE)
    }
    gam_daily <- function(df) {
      rhs <- paste(
        sprintf("s(tmax, k=%d)", cfg$k_tmax_day),
        sprintf("s(median_age, k=%d)", cfg$k_age),
        sprintf("s(humidity, k=%d)", cfg$k_hum),
        "median_income", "mean_college", "mean_hs",
        "mean_white", "mean_black", "mean_hispanic", "mean_asian",
        "mean_employed", "mean_unemployed", "mean_male", "mean_female",
        sprintf("s(date_num, bs='cs', k=%d)", cfg$k_time_date),
        sep = " + "
      )
      mgcv::gam(as.formula(paste0("CVD_per_100k ~ ", rhs)),
                data=df, family = mgcv::nb(link="log"), method="REML", select=TRUE)
    }
    gam_daily_ma <- function(df) {
      rhs <- paste(
        sprintf("s(tmax_3day_avg, k=%d)", cfg$k_tmax_day),
        sprintf("s(median_age, k=%d)", cfg$k_age),
        sprintf("s(humidity, k=%d)", cfg$k_hum),
        "median_income", "mean_college", "mean_hs",
        "mean_white", "mean_black", "mean_hispanic", "mean_asian",
        "mean_employed", "mean_unemployed", "mean_male", "mean_female",
        sprintf("s(date_num, bs='cs', k=%d)", cfg$k_time_date),
        sep = " + "
      )
      mgcv::gam(as.formula(paste0("CVD_per_100k ~ ", rhs)),
                data=df, family = mgcv::nb(link="log"), method="REML", select=TRUE)
    }
    
    pb$tick(sprintf("[%s] Fit: annual/monthly/daily/dailyMA", short))
    m_year  <- gam_yearly(py)
    m_month <- gam_monthly(pm)
    m_day   <- gam_daily(pd)
    m_dayma <- gam_daily_ma(pdma)
    
    # --- SAVE MODELS (RDS) ---
    rds_year  <- file.path(out_mod, sprintf("%s_%s_yearly.rds",  prefix, short))
    rds_month <- file.path(out_mod, sprintf("%s_%s_monthly.rds", prefix, short))
    rds_day   <- file.path(out_mod, sprintf("%s_%s_daily.rds",   prefix, short))
    rds_dayma <- file.path(out_mod, sprintf("%s_%s_dailyMA.rds", prefix, short))
    
    save_gam_with_means(m_year,  rds_year)
    save_gam_with_means(m_month, rds_month)
    save_gam_with_means(m_day,   rds_day)
    save_gam_with_means(m_dayma, rds_dayma)
    
    # --- COMPACT SUMMARY ROWS (one per scale) ---
    sum_year  <- cbind(outcome = short, scale = "yearly",  summarize_gam(m_year,  "mean_tmax"))
    sum_month <- cbind(outcome = short, scale = "monthly", summarize_gam(m_month, "mean_tmax"))
    sum_day   <- cbind(outcome = short, scale = "daily",   summarize_gam(m_day,   "^s\\(tmax"))
    sum_dayma <- cbind(outcome = short, scale = "dailyMA", summarize_gam(m_dayma, "tmax_3day_avg"))
    
    # --- WRITE/APPEND the summary table for the whole run ---
    summ_path <- file.path(out_tab, sprintf("gam_summary_%s_edf_p_r2_dev.csv", prefix))
    append_csv_row(sum_year,  summ_path)
    append_csv_row(sum_month, summ_path)
    append_csv_row(sum_day,   summ_path)
    append_csv_row(sum_dayma, summ_path)
    
    sink(file.path(out_log, paste0("gam_summaries_", outcome_type, "_", short, ".txt")))
    cat("\n==== ", toupper(outcome_type), " • ", short, " ====\n")
    print(summary(m_year));  print(gam.check(m_year))
    print(summary(m_month)); print(gam.check(m_month))
    print(summary(m_day));   print(gam.check(m_day))
    print(summary(m_dayma)); print(gam.check(m_dayma))
    sink()
    
    # Predict
    py$pred   <- predict(m_year,  type="response")
    pm$pred   <- predict(m_month, type="response")
    pd$pred   <- predict(m_day,   type="response")
    pdma$pred <- predict(m_dayma, type="response")
    
    pb$tick(sprintf("[%s] Scan thresholds (obs/pred; annual, monthly, daily, 3-day MA)", short))
    scan_obs_year  <- scan_threshold(py,   "mean_tmax", "CVD_per_100k", thr_seq = cfg$thr_year)
    scan_pred_year <- scan_threshold(py,   "mean_tmax", "CVD_per_100k", thr_seq = cfg$thr_year,  pred_vec = py$pred)
    scan_obs_mon   <- scan_threshold(pm,   "mean_tmax", "CVD_per_100k", thr_seq = cfg$thr_month)
    scan_pred_mon  <- scan_threshold(pm,   "mean_tmax", "CVD_per_100k", thr_seq = cfg$thr_month, pred_vec = pm$pred)
    scan_obs_day   <- scan_threshold(pd,   "tmax",      "CVD_per_100k", thr_seq = cfg$thr_day)
    scan_pred_day  <- scan_threshold(pd,   "tmax",      "CVD_per_100k", thr_seq = cfg$thr_day,   pred_vec = pd$pred)
    scan_obs_dayma <- scan_threshold(pdma, "tmax_3day_avg","CVD_per_100k", thr_seq = cfg$thr_day)
    scan_pred_dayma<- scan_threshold(pdma, "tmax_3day_avg","CVD_per_100k", thr_seq = cfg$thr_day, pred_vec = pdma$pred)
    
    # Save grids
    readr::write_csv(scan_obs_year$grid,   file.path(out_tab, paste0(prefix,"_",short,"_T1_threshold_grid_annual_observed.csv")))
    readr::write_csv(scan_pred_year$grid,  file.path(out_tab, paste0(prefix,"_",short,"_T2_threshold_grid_annual_predicted.csv")))
    readr::write_csv(scan_obs_mon$grid,    file.path(out_tab, paste0(prefix,"_",short,"_T3_threshold_grid_monthly_observed.csv")))
    readr::write_csv(scan_pred_mon$grid,   file.path(out_tab, paste0(prefix,"_",short,"_T4_threshold_grid_monthly_predicted.csv")))
    readr::write_csv(scan_obs_day$grid,    file.path(out_tab, paste0(prefix,"_",short,"_T5_threshold_grid_daily_observed.csv")))
    readr::write_csv(scan_pred_day$grid,   file.path(out_tab, paste0(prefix,"_",short,"_T6_threshold_grid_daily_predicted.csv")))
    readr::write_csv(scan_obs_dayma$grid,  file.path(out_tab, paste0(prefix,"_",short,"_T7_threshold_grid_dailyMA_observed.csv")))
    readr::write_csv(scan_pred_dayma$grid, file.path(out_tab, paste0(prefix,"_",short,"_T8_threshold_grid_dailyMA_predicted.csv")))
    
    # Peaks
    collect_peaks <- function(nm, scan) tibble::tibble(
      scale = nm,
      peak_threshold = scan$peak$threshold,
      peak_avg_excess = scan$peak$avg_excess_rate,
      highest_pos_threshold = scan$highest_positive$threshold,
      highest_pos_avg_excess = scan$highest_positive$avg_excess_rate
    )
    peaks <- dplyr::bind_rows(
      collect_peaks("annual_observed",  scan_obs_year),
      collect_peaks("annual_predicted", scan_pred_year),
      collect_peaks("monthly_observed", scan_obs_mon),
      collect_peaks("monthly_predicted",scan_pred_mon),
      collect_peaks("daily_observed",   scan_obs_day),
      collect_peaks("daily_predicted",  scan_pred_day),
      collect_peaks("dailyMA_observed", scan_obs_dayma),
      collect_peaks("dailyMA_predicted",scan_pred_dayma)
    )
    readr::write_csv(peaks, file.path(out_tab, paste0(prefix,"_",short,"_S0_threshold_peaks.csv")))
    
    # ---- Bootstrap CI bands (progress inside each call) ----
    pb$tick(sprintf("[%s] Bootstrap CIs (annual obs/pred)", short))
    ci_ann_obs <- bootstrap_threshold_curve(py,   temp_var = "mean_tmax", pred_var = "CVD_per_100k", thr_seq = cfg$thr_year,  R = cfg$boot_R)
    ci_ann_pre <- bootstrap_threshold_curve(py,   temp_var = "mean_tmax", pred_var = "pred",          thr_seq = cfg$thr_year,  R = cfg$boot_R)
    
    pb$tick(sprintf("[%s] Bootstrap CIs (monthly obs/pred)", short))
    ci_mon_obs <- bootstrap_threshold_curve(pm,   temp_var = "mean_tmax", pred_var = "CVD_per_100k", thr_seq = cfg$thr_month, R = cfg$boot_R)
    ci_mon_pre <- bootstrap_threshold_curve(pm,   temp_var = "mean_tmax", pred_var = "pred",          thr_seq = cfg$thr_month, R = cfg$boot_R)
    
    pb$tick(sprintf("[%s] Bootstrap CIs (daily obs/pred)", short))
    ci_day_obs <- bootstrap_threshold_curve(pd,   temp_var = "tmax",      pred_var = "CVD_per_100k", thr_seq = cfg$thr_day,   R = cfg$boot_R)
    ci_day_pre <- bootstrap_threshold_curve(pd,   temp_var = "tmax",      pred_var = "pred",          thr_seq = cfg$thr_day,   R = cfg$boot_R)
    
    pb$tick(sprintf("[%s] Bootstrap CIs (3-day MA obs/pred)", short))
    ci_dma_obs <- bootstrap_threshold_curve(pdma, temp_var = "tmax_3day_avg", pred_var = "CVD_per_100k", thr_seq = cfg$thr_day, R = cfg$boot_R)
    ci_dma_pre <- bootstrap_threshold_curve(pdma, temp_var = "tmax_3day_avg", pred_var = "pred",          thr_seq = cfg$thr_day, R = cfg$boot_R)
    
    # Merge CIs to grids
    g_ann_obs <- dplyr::left_join(scan_obs_year$grid,  ci_ann_obs, by = "threshold")
    g_ann_pre <- dplyr::left_join(scan_pred_year$grid, ci_ann_pre, by = "threshold")
    g_mon_obs <- dplyr::left_join(scan_obs_mon$grid,   ci_mon_obs, by = "threshold")
    g_mon_pre <- dplyr::left_join(scan_pred_mon$grid,  ci_mon_pre, by = "threshold")
    g_day_obs <- dplyr::left_join(scan_obs_day$grid,   ci_day_obs, by = "threshold")
    g_day_pre <- dplyr::left_join(scan_pred_day$grid,  ci_day_pre, by = "threshold")
    g_dma_obs <- dplyr::left_join(scan_obs_dayma$grid, ci_dma_obs, by = "threshold")
    g_dma_pre <- dplyr::left_join(scan_pred_dayma$grid,ci_dma_pre, by = "threshold")
    
    # Plots
    pb$tick(sprintf("[%s] Plotting (8 figs)", short))
    plot_thr(g_ann_obs, scan_obs_year$peak,  scan_obs_year$highest_positive,
             "Mean Tmax (°C), annual (observed)",
             file.path(out_fig, paste0(prefix,"_",short,"_obs_annual_excess.png")),
             ci_lo="lower_ci", ci_hi="upper_ci")
    plot_thr(g_ann_pre, scan_pred_year$peak, scan_pred_year$highest_positive,
             "Mean Tmax (°C), annual (predicted)",
             file.path(out_fig, paste0(prefix,"_",short,"_pred_annual_excess.png")),
             ci_lo="lower_ci", ci_hi="upper_ci")
    plot_thr(g_mon_obs, scan_obs_mon$peak,   scan_obs_mon$highest_positive,
             "Mean Tmax (°C), monthly (observed)",
             file.path(out_fig, paste0(prefix,"_",short,"_obs_monthly_excess.png")),
             ci_lo="lower_ci", ci_hi="upper_ci")
    plot_thr(g_mon_pre, scan_pred_mon$peak,  scan_pred_mon$highest_positive,
             "Mean Tmax (°C), monthly (predicted)",
             file.path(out_fig, paste0(prefix,"_",short,"_pred_monthly_excess.png")),
             ci_lo="lower_ci", ci_hi="upper_ci")
    plot_thr(g_day_obs, scan_obs_day$peak,   scan_obs_day$highest_positive,
             "Daily Tmax (°C), observed",
             file.path(out_fig, paste0(prefix,"_",short,"_obs_daily_excess.png")),
             ci_lo="lower_ci", ci_hi="upper_ci")
    plot_thr(g_day_pre, scan_pred_day$peak,  scan_pred_day$highest_positive,
             "Daily Tmax (°C), predicted",
             file.path(out_fig, paste0(prefix,"_",short,"_pred_daily_excess.png")),
             ci_lo="lower_ci", ci_hi="upper_ci")
    plot_thr(g_dma_obs, scan_obs_dayma$peak, scan_obs_dayma$highest_positive,
             "3-day mean Tmax (°C), observed",
             file.path(out_fig, paste0(prefix,"_",short,"_obs_dailyMA_excess.png")),
             ci_lo="lower_ci", ci_hi="upper_ci")
    plot_thr(g_dma_pre, scan_pred_dayma$peak,scan_pred_dayma$highest_positive,
             "3-day mean Tmax (°C), predicted",
             file.path(out_fig, paste0(prefix,"_",short,"_pred_dailyMA_excess.png")),
             ci_lo="lower_ci", ci_hi="upper_ci")
    
    # Map
    pb$tick(sprintf("[%s] Mapping peak annual threshold", short))
    thr_ann <- (peaks %>% dplyr::filter(scale=="annual_predicted"))$peak_threshold[1]
    if (is.finite(thr_ann)) {
      yr_map <- py %>%
        dplyr::group_by(community) %>%
        dplyr::mutate(heat_year = mean_tmax > thr_ann) %>%
        dplyr::summarise(
          mean_pred_hot  = mean(pred[heat_year],  na.rm=TRUE),
          mean_pred_cool = mean(pred[!heat_year], na.rm=TRUE),
          excess_rate_per_year = mean_pred_hot - mean_pred_cool,
          n_heat_years = sum(heat_year),
          total_excess = excess_rate_per_year * n_heat_years,
          .groups="drop"
        )
      map_dat <- dplyr::left_join(cas, yr_map, by="community")
      
      x <- map_dat$excess_rate_per_year
      scale_info <- safe_diverging_scale(x)
      if (!scale_info$ok) {
        message(sprintf("  ↳ [%s] Map skipped (no finite excess values).", short))
      } else {
        png(file.path(out_fig, paste0(prefix,"_",short,"_pred_annual_rate_map.png")),
            width=1400, height=1000, res=200)
        print(
          tm_shape(map_dat) +
            tm_polygons(
              "excess_rate_per_year",
              fill.scale = tm_scale_continuous(
                values = "-brewer.RdBu",     # blue for −, red for +
                midpoint = 0,
                limits   = scale_info$limits,
                ticks    = scale_info$ticks
              ),
              fill.legend = tm_legend(
                title = paste0("Excess ", short, " rate per 100k\n(annual; peak thr)"),
                reverse = FALSE
              ),
              col = "white", col_alpha = 0.5, lwd = 0.4
            ) +
            tm_layout(frame = FALSE, legend.outside = TRUE, legend.outside.position = "right")
        )
        dev.off()
      }
    }
    
    # Sensitivity (≤2019)
    pb$tick(sprintf("[%s] Sensitivity (≤2019)", short))
    if (cfg$Covid_exclusion) {
      yr19 <- py %>% dplyr::filter(year <= 2019)
      if (nrow(yr19)) {
        m_year_19 <- mgcv::gam(update(formula(m_year), . ~ .), data = yr19,
                               family = mgcv::nb(link="log"), method="REML", select=TRUE)
        yr19$pred <- predict(m_year_19, type="response")
        scan_19 <- scan_threshold(yr19, "mean_tmax", "CVD_per_100k", thr_seq = cfg$thr_year, pred_vec = yr19$pred)
        readr::write_csv(scan_19$grid, file.path(out_tab, paste0(prefix,"_",short,"_S1_thresholds_annual_covid_excl.csv")))
      }
    }
  } # end endpoint loop
  
  pb$end()
  invisible(TRUE)
}



# -----------------------
# Run both outcomes
# -----------------------
run_pipeline("mortality")
run_pipeline("ed")

# -----------------------
# Plot the relationship using the model output
# -----------------------

prefix <- "ED"

# optional, but great for avoiding label overlap:
suppressPackageStartupMessages({ library(ggrepel) })

summs <- readr::read_csv(file.path(out_tab, sprintf("gam_summary_%s_edf_p_r2_dev.csv", prefix)),
                         show_col_types = FALSE) %>%
  mutate(
    outcome = factor(outcome, levels = c("CVD","CHD","MI","Stroke")),
    scale_raw = factor(scale, levels = c("yearly","monthly","daily","dailyMA")),
    scale = recode(scale_raw,
                   "yearly"  = "Warm Months",
                   "monthly" = "Month",
                   "daily"   = "Day",
                   "dailyMA" = "Daily (lag 0–3)"),
    scale = factor(scale, levels = c("Warm Months","Month","Day","Daily (lag 0–3)")),
    p_label = dplyr::case_when(
      is.na(p)    ~ "p=NA",
      p < 0.001   ~ "p<0.001",
      TRUE        ~ paste0("p=", format(round(p, 3), nsmall = 3))
    )
  )

pd <- position_dodge(width = 0.45)

p_bw <- ggplot(summs, aes(x = outcome, y = edf, shape = scale)) +
  geom_hline(yintercept = 1, linetype = 3, linewidth = 0.4, color = "grey40") +
  geom_point(size = 3.3, position = pd, color = "black") +
  ggrepel::geom_text_repel(aes(label = p_label),
                           position = pd, size = 3.2, min.segment.length = 0.05,
                           box.padding = 0.2, point.padding = 0.2, seed = 123) +
  scale_shape_manual(values = c(16, 17, 15, 3)) +  # circle, triangle, square, plus
  labs(
    x = NULL,
    y = "Nonlinearity (edf of s(Tmax))",
    title = "Heat–mortality smooth: nonlinearity & significance",
    subtitle = "Shapes denote time scale; dashed line marks ~linear (edf≈1)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(file.path(out_fig, sprintf("%s_compare_sTmax_edf_p_bw.png", prefix)),
       p_bw, width = 7.2, height = 4.6, dpi = 600)
ggsave(file.path(out_fig, sprintf("%s_compare_sTmax_edf_p_bw.pdf", prefix)),
       p_bw, width = 7.2, height = 4.6, device = cairo_pdf)




# read the GAM edf/p table
summs <- read_csv(file.path(out_tab, sprintf("gam_summary_%s_edf_p_r2_dev.csv", prefix)),
                  show_col_types = FALSE) %>%
  mutate(
    outcome = factor(outcome, levels = c("CVD","CHD","MI","Stroke")),
    scale   = recode(scale,
                     "yearly"  = "Warm Months",
                     "monthly" = "Month",
                     "daily"   = "Day",
                     "dailyMA" = "Daily (lag 0–3)"),
    scale   = factor(scale, levels = c("Warm Months","Month","Day","Daily (lag 0–3)")),
    mlogp   = -log10(p)
  )

summs <- summs %>%
  mutate(sig = cut(p,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("p<0.001","p<0.01","p<0.05","NS"),
                   right = TRUE))

p2 <- ggplot(summs, aes(x = scale, y = edf, fill = sig)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_col(width = 0.6, color = "grey25") +
  facet_wrap(~ outcome, nrow = 1, scales = "fixed") +
  scale_fill_manual(
    values = c(
      "p<0.001" = "#B2182B",
      "p<0.01"  = "#D6604D",
      "p<0.05"  = "#F4A582",
      "NS"      = "#6C757D"
    ),
    name = "Significance"
  ) +
  labs(
    x = NULL,
    y = "Curvature of the temperature–mortality\nrelationship (edf)",
    title = "Strength and Shape of the Temperature–Morbidity\nRelationship Across Time Scales"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p2

ggsave(file.path(out_fig, sprintf("barchart_compare_sTmax_edf_p_color.png", prefix)),
       p2, width = 9.2, height = 8.6, dpi = 600)

# make the plot
p_tile <- ggplot(summs, aes(scale, outcome, fill = mlogp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("edf = %.2f", edf)), size = 3.1) +
  scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey90",
                      name = expression(-log[10](p))) +
  labs(
    x = NULL, y = NULL,
    title = "Significance (–log10 p) of heat–mortality smooths",
    subtitle = "edf (effective degrees of freedom) annotated within each cell"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

# save high-resolution PNG (change to PDF/TIFF if desired)
ggsave(
  filename = file.path(out_fig, sprintf("%s_tile_sTmax_significance.png", prefix)),
  plot = p_tile,
  width = 7.5, height = 4.5, dpi = 600
)


# -----------------------
# Tabulating the peak threshold for CVD and its subtypes by outcome type
# -----------------------

# ---- PATH ----
path_tables <- "/Users/saborpete/Desktop/Peter/PhD/Dissertation/Analysis/Code and Data/Heat_CVD_Chicago/output/tables"

files <- list.files(path_tables, full.names = TRUE, pattern = "\\.csv$")

# ---- Helper to extract subtype + temporal scale from filename ----
parse_metadata <- function(fname) {
  
  base <- basename(fname)
  
  tibble(
    source_file = base,
    # ED vs MORT
    type = case_when(
      str_detect(base, "^ED_") ~ "ED",
      str_detect(base, "^MORT") ~ "MORT",
      TRUE ~ NA_character_
    ),
    # subtype is the piece after ED_ or MORT_ until next underscore
    subtype = base %>%
      str_remove("^ED_") %>%
      str_remove("^MORT_") %>%
      str_extract("^[A-Za-z0-9]+"),
    # predicted?
    predicted_flag = str_detect(base, "predicted"),
    # temporal scale
    temporal = case_when(
      str_detect(base, "annual")  ~ "annual",
      str_detect(base, "monthly") ~ "monthly",
      str_detect(base, "dailyMA") ~ "dailyMA",
      str_detect(base, "daily")   ~ "daily",
      TRUE ~ NA_character_
    )
  )
}

# ---- Load with metadata ----
load_file <- function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  meta <- parse_metadata(f)
  
  df %>%
    mutate(
      source_file = meta$source_file,
      type        = meta$type,
      subtype     = meta$subtype,
      temporal    = meta$temporal,
      predicted   = meta$predicted_flag
    )
}

dat <- map_dfr(files, load_file)

# ---- Keep only predicted ----
dat_pred <- dat %>% filter(predicted)

# ---- Summary function ----
compute_summary <- function(df) {
  
  peak <- df %>%
    filter(avg_excess_rate == max(avg_excess_rate, na.rm = TRUE)) %>%
    slice(1)
  
  zero <- df %>%
    arrange(threshold) %>%
    filter(avg_excess_rate <= 0) %>%
    slice(1)
  
  tibble(
    type        = df$type[1],
    subtype     = df$subtype[1],
    temporal    = df$temporal[1],
    
    peak_threshold       = peak$threshold,
    peak_excess_rate     = peak$avg_excess_rate,
    peak_sd              = peak$sd_excess_rate,
    
    zero_cross_threshold = zero$threshold,
    zero_cross_excess    = zero$avg_excess_rate,
    zero_cross_sd        = zero$sd_excess_rate
  )
}

# ---- Build full results table ----
results <- dat_pred %>%
  group_by(type, subtype, temporal) %>%
  group_modify(~ compute_summary(.x)) %>%
  ungroup() %>%
  arrange(type, subtype, temporal)

# ---- WRITE OUT FULL TABLE ----
write_csv(results, file.path(path_tables, "results_all.csv"))

# ---- SPLIT AND WRITE OUT ----
results_ED <- results %>% filter(type == "ED")
results_MORT <- results %>% filter(type == "MORT")

write_csv(results_ED,   file.path(path_tables, "results_ED.csv"))
write_csv(results_MORT, file.path(path_tables, "results_MORT.csv"))


# -----------------------
# Plotting the peak thresholds and excess rates
# -----------------------


# assuming `results` is already in your environment
results_plot <- results %>%
  mutate(
    type = factor(type,
                  levels = c("MORT", "ED"),
                  labels = c("Mortality", "ED visits")),
    temporal = factor(
      temporal,
      levels = c("annual", "monthly", "daily", "dailyMA"),
      labels = c("Annual", "Monthly", "Daily", "Daily (MA)")
    ),
    # tweak order/labels as you like
    subtype = factor(
      subtype,
      levels = c("CVD", "CHD", "MI", "Stroke"),
      labels = c("All CVD", "CHD", "MI", "Stroke")
    )
  )

# long format for thresholds
thresholds_long <- results_plot %>%
  select(type, subtype, temporal,
         peak_threshold, zero_cross_threshold) %>%
  pivot_longer(
    cols = c(peak_threshold, zero_cross_threshold),
    names_to = "kind",
    values_to = "threshold"
  ) %>%
  mutate(
    kind = recode(kind,
                  peak_threshold       = "Peak excess",
                  zero_cross_threshold = "Zero-crossing")
  )



pd_thr <- position_dodge(width = 0.45)

p_thresh <- ggplot(thresholds_long,
                   aes(x = threshold, y = subtype,
                       color = temporal, shape = kind)) +
  # reference "comfortable" temperature if you want (e.g. 21°C)
  geom_vline(xintercept = 21,
             linetype = 2, linewidth = 0.4, color = "grey40") +
  geom_point(size = 3.1, position = pd_thr) +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  scale_color_brewer(palette = "Set1",
                     name = "Temporal scale") +
  scale_shape_manual(
    name   = "Threshold type",
    values = c(`Peak excess` = 16, `Zero-crossing` = 17)
  ) +
  labs(
    x = "Temperature threshold (°C)",
    y = NULL,
    title = "Heat thresholds for CVD endpoints in Chicago",
    subtitle = "Peak excess and zero-crossing thresholds by temporal scale and outcome type"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    legend.box         = "vertical",
    strip.text         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 11),
    axis.text.y        = element_text(size = 10)
  )

p_thresh

# save
ggsave(
  file.path(path_tables, "fig_thresh_MORT_ED_cvd_subtypes.png"),
  p_thresh, width = 7.2, height = 6.0, dpi = 600
)


pd_peak <- position_dodge(width = 0.35)

p_peak <- ggplot(results_plot,
                 aes(x = temporal, y = peak_excess_rate,
                     color = type, shape = type, group = type)) +
  geom_hline(yintercept = 0, linetype = 3,
             linewidth = 0.4, color = "grey40") +
  geom_errorbar(
    aes(ymin = peak_excess_rate - peak_sd,
        ymax = peak_excess_rate + peak_sd),
    position = pd_peak, width = 0,
    linewidth = 0.6
  ) +
  geom_point(position = pd_peak, size = 2.8) +
  facet_wrap(~ subtype, nrow = 1) +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  scale_shape_manual(
    name   = NULL,
    values = c(`Mortality` = 16, `ED visits` = 17)
  ) +
  labs(
    x = NULL,
    y = "Peak excess rate (per 100,000, or units used)",
    title = "Magnitude of peak excess rates across temporal scales",
    subtitle = "Comparison of mortality vs ED visits by CVD subtype"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    legend.box         = "horizontal",
    strip.text         = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 25, hjust = 1),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 11)
  )

p_peak

ggsave(
  file.path(path_tables, "fig_peak_excess_MORT_ED_cvd_subtypes.png"),
  p_peak, width = 7.5, height = 3.8, dpi = 600
)


# Peak thresholds only
peak_thr <- results_plot %>%
  select(type, subtype, temporal, peak_threshold)

p_thresh_clean <- ggplot(peak_thr,
                         aes(x = peak_threshold, y = subtype, color = subtype)) +
  geom_bar() +
  facet_grid(type ~ temporal) +
  labs(
    x = "Peak temperature threshold (°C)",
    y = NULL,
    title = "Heat thresholds for CVD endpoints in Chicago",
    subtitle = "Peak excess thresholds by temporal scale, outcome type, and CVD subtype"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    strip.text         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 11),
    axis.text.y        = element_text(size = 10)
  )

p_thresh_clean

ggsave(
  file.path(path_tables, "fig_thresh_peak_only_MORT_ED_cvd_subtypes.png"),
  p_thresh_clean, width = 7.2, height = 5.5, dpi = 600
)


# Compute a sensible y-range based on central 90% of (peak ± SD)
peak_range <- results_plot %>%
  transmute(
    lo = peak_excess_rate - peak_sd,
    hi = peak_excess_rate + peak_sd
  ) %>%
  summarise(
    lo_q = quantile(lo, 0.05, na.rm = TRUE),
    hi_q = quantile(hi, 0.95, na.rm = TRUE)
  )

y_lo <- as.numeric(peak_range$lo_q)
y_hi <- as.numeric(peak_range$hi_q)

pd_peak <- position_dodge(width = 0.35)

p_peak_clean <- ggplot(results_plot,
                       aes(x = temporal, y = peak_excess_rate,
                           color = type, shape = type, group = type)) +
  geom_hline(yintercept = 0, linetype = 3,
             linewidth = 0.4, color = "grey40") +
  geom_errorbar(
    aes(ymin = peak_excess_rate - peak_sd,
        ymax = peak_excess_rate + peak_sd),
    position = pd_peak, width = 0,
    linewidth = 0.6
  ) +
  geom_point(position = pd_peak, size = 2.8) +
  facet_wrap(~ subtype, ncol = 2) +
  coord_cartesian(ylim = c(y_lo, y_hi)) +
  scale_color_manual(
    name   = NULL,
    values = c("Mortality" = "#1b9e77", "ED visits" = "#d95f02")
  ) +
  scale_shape_manual(
    name   = NULL,
    values = c("Mortality" = 16, "ED visits" = 17)
  ) +
  labs(
    x = NULL,
    y = "Peak excess rate (per 100,000, or units used)",
    title = "Magnitude of peak excess rates across temporal scales",
    subtitle = "Comparison of mortality vs ED visits by CVD subtype"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    legend.box         = "horizontal",
    strip.text         = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 25, hjust = 1),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 11)
  )

p_peak_clean

ggsave(
  file.path(path_tables, "fig_peak_excess_zoom_MORT_ED_cvd_subtypes.png"),
  p_peak_clean, width = 7.0, height = 5.0, dpi = 600
)


peak_thr <- results_plot %>%
  select(type, subtype, temporal, peak_threshold)

p_thresh_bar <- ggplot(peak_thr,
                       aes(x = subtype, y = peak_threshold, fill = subtype)) +
  geom_col(width = 0.65) +
  facet_grid(type ~ temporal) +
  labs(
    x = NULL,
    y = "Peak temperature threshold (°C)",
    title = "Heat thresholds for CVD endpoints in Chicago",
    subtitle = "Peak excess thresholds by temporal scale, outcome type, and CVD subtype"
  ) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11)
  )

p_thresh_bar

ggsave(
  file.path(path_tables, "fig_thresh_peak_bar_MORT_ED_cvd_subtypes.png"),
  p_thresh_bar, width = 8, height = 5.5, dpi = 600
)


# Rename temporal scale levels
results_plot2 <- results_plot %>%
  mutate(
    temporal = fct_recode(
      temporal,
      "Warm Months"     = "Annual",
      "Daily (lag 0–3)" = "Daily (MA)"
    )
  )

# Recompute y-range
peak_range <- results_plot2 %>%
  transmute(
    lo = peak_excess_rate - peak_sd,
    hi = peak_excess_rate + peak_sd
  ) %>%
  summarise(
    lo_q = quantile(lo, 0.05, na.rm = TRUE),
    hi_q = quantile(hi, 0.95, na.rm = TRUE)
  )

y_lo <- as.numeric(peak_range$lo_q)
y_hi <- as.numeric(peak_range$hi_q)

pd_peak <- position_dodge(width = 0.35)

p_peak_clean2 <- ggplot(results_plot2,
                        aes(x = temporal, y = peak_excess_rate,
                            color = type, shape = type, group = type)) +
  geom_hline(yintercept = 0, linetype = 3,
             linewidth = 0.4, color = "grey40") +
  geom_errorbar(
    aes(ymin = peak_excess_rate - peak_sd,
        ymax = peak_excess_rate + peak_sd),
    position = pd_peak, width = 0, linewidth = 0.6
  ) +
  geom_point(position = pd_peak, size = 2.8) +
  facet_wrap(~ subtype, ncol = 2) +
  coord_cartesian(ylim = c(y_lo, y_hi)) +
  scale_color_manual(
    name = NULL,
    values = c("Mortality" = "#1b9e77", "ED visits" = "#d95f02")
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Mortality" = 16, "ED visits" = 17)
  ) +
  labs(
    x = NULL,
    y = "Peak excess rate (per 100,000 population)",
    title = "Magnitude of peak excess rates across temporal scales",
    subtitle = "Comparison of mortality vs emergency department visits by cardiovascular disease subtype"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

p_peak_clean2

ggsave(
  file.path(path_tables, "fig_peak_excess_zoom_renamed_MORT_ED_cvd_subtypes.png"),
  p_peak_clean2, width = 7.3, height = 5.0, dpi = 600
)



# Rename temporal scale levels using Set A
peak_thr2 <- results_plot %>%
  mutate(
    temporal = fct_recode(
      temporal,
      "Warm Months"     = "Annual",
      "Monthly"         = "Monthly",
      "Daily"           = "Daily",
      "Daily (lag 0–3)" = "Daily (MA)"
    )
  ) %>%
  select(type, subtype, temporal, peak_threshold)

p_thresh_bar2 <- ggplot(peak_thr2,
                        aes(x = subtype, y = peak_threshold, fill = subtype)) +
  geom_col(width = 0.65) +
  facet_grid(type ~ temporal) +
  labs(
    x = NULL,
    y = "Peak temperature threshold (°C)",
    title = "Heat thresholds for CVD endpoints in Chicago",
    subtitle = "Peak excess thresholds by temporal scale, outcome type, and CVD subtype"
  ) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11)
  )

p_thresh_bar2

ggsave(
  file.path(path_tables, "fig_thresh_peak_bar_renamed_MORT_ED_cvd_subtypes.png"),
  p_thresh_bar2, width = 8, height = 5.5, dpi = 600
)


## COVID Sensitivity

# --- metadata parser: type, subtype, temporal, scenario ---
parse_meta_sens <- function(fname) {
  base <- basename(fname)
  
  tibble(
    source_file = base,
    type = case_when(
      str_detect(base, "^ED_")   ~ "ED",
      str_detect(base, "^MORT") ~ "MORT",
      TRUE ~ NA_character_
    ),
    subtype = base %>%
      str_remove("^ED_") %>%
      str_remove("^MORT_") %>%
      str_extract("^[A-Za-z0-9]+"),
    temporal = case_when(
      str_detect(base, "annual")  ~ "annual",
      str_detect(base, "monthly") ~ "monthly",
      str_detect(base, "dailyMA") ~ "dailyMA",
      str_detect(base, "daily")   ~ "daily",
      TRUE ~ NA_character_
    ),
    scenario = case_when(
      str_detect(base, "_S1_") | str_detect(base, "covid_excl") ~ "covid_excluded",
      TRUE ~ "main"
    ),
    # keep main predicted + covid_excl files
    predicted_flag = str_detect(base, "predicted") |
      str_detect(base, "covid_excl")
  )
}

load_file_sens <- function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  meta <- parse_meta_sens(f)
  
  df %>%
    mutate(
      source_file = meta$source_file,
      type        = meta$type,
      subtype     = meta$subtype,
      temporal    = meta$temporal,
      scenario    = meta$scenario,
      predicted   = meta$predicted_flag
    )
}

dat_all <- map_dfr(files, load_file_sens)

# keep only predicted / covid-excluded thresholds
dat_pred <- dat_all %>%
  filter(predicted, !is.na(type), !is.na(subtype), !is.na(temporal))

# --- function to compute peak threshold for each group ---
compute_summary_sens <- function(df) {
  peak <- df %>%
    filter(avg_excess_rate == max(avg_excess_rate, na.rm = TRUE)) %>%
    slice(1)
  
  tibble(
    type           = df$type[1],
    subtype        = df$subtype[1],
    temporal       = df$temporal[1],
    scenario       = df$scenario[1],
    peak_threshold = peak$threshold,
    peak_excess    = peak$avg_excess_rate
  )
}

results_sens <- dat_pred %>%
  group_by(type, subtype, temporal, scenario) %>%
  group_modify(~ compute_summary_sens(.x)) %>%
  ungroup()

covid_peaks <- results_sens %>%
  filter(temporal == "annual") %>%      # Warm Months
  mutate(
    type = recode(type,
                  "MORT" = "Mortality",
                  "ED"   = "ED visits"),
    subtype = factor(
      subtype,
      levels = c("CVD", "CHD", "MI", "Stroke"),
      labels = c("All CVD", "CHD", "MI", "Stroke")
    ),
    scenario = factor(
      scenario,
      levels = c("main", "covid_excluded"),
      labels = c("Main analysis", "COVID period excluded")
    )
  )

p_covid <- ggplot(covid_peaks,
                  aes(x = subtype, y = peak_threshold, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.55) +
  facet_wrap(~ type, nrow = 1) +
  labs(
    x = NULL,
    y = "Peak temperature threshold (°C)",
    title = "Warm-season heat thresholds in COVID sensitivity analysis",
    subtitle = "Peak excess thresholds with and without COVID period, by CVD subtype and outcome type",
    fill = NULL
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 25, hjust = 1),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 11),
    legend.position  = "top"
  )

p_covid

ggsave(
  file.path(path_tables, "fig_covid_sensitivity_warmmonths_thresholds.png"),
  p_covid, width = 7.3, height = 4.6, dpi = 600
)


####### ############################################# PCA + k-means for community area clustering

# -----------------------------
# 1) Build mortality panel_yearly_simple
# -----------------------------

# Read deaths and add DATE
deaths_raw <- try_read_csv(paths$deaths)
stopifnot(!is.null(deaths_raw))

# Derive DATE column (DOD or DATE)
date_col <- c("DOD","DATE")[c("DOD","DATE") %in% names(deaths_raw)][1]
stopifnot(!is.null(date_col))

deaths <- deaths_raw %>%
  dplyr::mutate(DATE = as.Date(.data[[date_col]]))

# Make sure we have 'community'
deaths <- ensure_community(deaths)
if (!"community" %in% names(deaths)) {
  stop("Outcome file: could not make 'community' column.")
}

# Which outcome flags are present?
present <- intersect(expected_flags, names(deaths))
if (!length(present)) stop("No CVD/CHD/MI/Stroke flags found in the mortality file.")

# Warm months + study years
deaths <- deaths %>%
  dplyr::filter(lubridate::month(DATE) %in% 5:9) %>%
  dplyr::filter(lubridate::year(DATE) >= cfg$start_year,
                lubridate::year(DATE) <= cfg$end_year)

# ---- Daily → yearly counts by community ----
counts_daily_simple <- deaths %>%
  dplyr::group_by(community, DATE) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(present), ~ sum(.x, na.rm = TRUE), .names = "{.col}_count"),
    all_deaths = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    year  = lubridate::year(DATE),
    month = lubridate::month(DATE),
    day   = lubridate::day(DATE)
  )

counts_yearly_simple <- counts_daily_simple %>%
  dplyr::group_by(community, year) %>%
  dplyr::summarise(
    dplyr::across(dplyr::ends_with("_count"), sum, na.rm = TRUE),
    all_deaths = sum(all_deaths, na.rm = TRUE),
    .groups = "drop"
  )

# ---- Join ACS + temps_yearly and compute rates ----
panel_yearly_simple <- acs_summary %>%
  dplyr::inner_join(counts_yearly_simple, by = c("community","year")) %>%
  dplyr::inner_join(temps_yearly,        by = c("community","year")) %>%
  mk_rates()

# Add NO2, PM2.5, NDVI (z-scored)
if (!is.null(pm25_tbl)) {
  panel_yearly_simple <- panel_yearly_simple %>%
    dplyr::left_join(pm25_tbl, by = "year")
}
if (!is.null(no2_tbl)) {
  panel_yearly_simple <- panel_yearly_simple %>%
    dplyr::left_join(no2_tbl, by = "year")
}
if (!is.null(ndvi_year)) {
  panel_yearly_simple <- panel_yearly_simple %>%
    dplyr::left_join(ndvi_year, by = c("community","year"))
}

panel_yearly_simple <- z_air_green(panel_yearly_simple)

# -----------------------------
# 2) Build community-level feature set for PCA
# -----------------------------

feat <- panel_yearly_simple %>%
  dplyr::filter(year >= 2011, year <= 2022) %>%
  dplyr::group_by(community) %>%
  dplyr::summarise(
    # Climate
    mean_tmax        = mean(mean_tmax,      na.rm = TRUE),
    mean_humidity    = mean(mean_humidity,  na.rm = TRUE),
    
    # Demographics / SES
    median_income    = mean(median_income,  na.rm = TRUE),
    median_age       = mean(median_age,     na.rm = TRUE),
    mean_college     = mean(mean_college,   na.rm = TRUE),
    mean_hs          = mean(mean_hs,        na.rm = TRUE),
    mean_white       = mean(mean_white,     na.rm = TRUE),
    mean_black       = mean(mean_black,     na.rm = TRUE),
    mean_hispanic    = mean(mean_hispanic,  na.rm = TRUE),
    mean_asian       = mean(mean_asian,     na.rm = TRUE),
    mean_employed    = mean(mean_employed,  na.rm = TRUE),
    mean_unemployed  = mean(mean_unemployed,na.rm = TRUE),
    
    # Long-term environmental (already z-scored by z_air_green)
    pm25             = mean(pm25,           na.rm = TRUE),
    no2              = mean(no2,            na.rm = TRUE),
    ndvi             = mean(ndvi,           na.rm = TRUE),
    
    .groups = "drop"
  )

# Matrix for PCA
feat_mat <- as.matrix(feat %>% 
                        select(-community, -no2, -pm25))  # remove constant-air-pollution variables

# Drop zero-variance columns generally
sd_vals <- apply(feat_mat, 2, sd, na.rm = TRUE)
feat_mat_filtered <- feat_mat[, sd_vals > 0]

feat_scaled <- scale(feat_mat_filtered)

# -----------------------------
# 3) PCA + k-means clustering
# -----------------------------

pca <- prcomp(feat_scaled, center = FALSE, scale. = FALSE)

# Variance explained (for methods text)
summary(pca)

# Keep first 3 PCs
pc_scores <- tibble::as_tibble(pca$x[, 1:3], rownames = "community") %>%
  dplyr::rename(PC1 = PC1, PC2 = PC2, PC3 = PC3)

# k-means clustering on PC scores
set.seed(2025)
K <- 3  # or 4, depending on what looks interpretable
km <- kmeans(pc_scores %>% dplyr::select(PC1, PC2, PC3),
             centers = K, nstart = 100)

clusters <- pc_scores %>%
  dplyr::mutate(cluster = factor(km$cluster))

# -----------------------------
# 4) Warm-season CVD rate (2011–2019) per community
# -----------------------------

cvd_by_comm <- panel_yearly_simple %>%
  dplyr::filter(year >= 2011, year <= 2022) %>%
  dplyr::group_by(community) %>%
  dplyr::summarise(
    mean_CVD_per_100k = mean(CVD_per_100k, na.rm = TRUE),
    .groups = "drop"
  )

clusters <- tibble::tibble(
  community = feat$community,         # character names (ALBANY PARK, etc.)
  cluster   = factor(km$cluster) # 1,2,3 → factor
)

cluster_excess <- clusters %>%
  dplyr::left_join(cvd_by_comm, by = "community")

loading_tbl <- pca$rotation %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable")

readr::write_csv(loading_tbl,
                 file.path(out_tab, "PCA_loadings_community_characteristics.csv"))

map_clusters <- cas %>%
  dplyr::left_join(clusters, by = "community")

tmap::tmap_mode("plot")

p_cluster_map <- tm_shape(map_clusters) +
  tm_polygons("cluster",
              palette = "Set1",
              title = "Community cluster") +
  tm_layout(frame = FALSE,
            legend.outside = TRUE,
            legend.outside.position = "right")

tmap::tmap_save(
  p_cluster_map,
  filename = file.path(out_fig, "cluster_map_community_chars.png"),
  width = 7, height = 6, dpi = 600
)

p_box <- ggplot(cluster_excess,
                aes(x = cluster, y = mean_CVD_per_100k)) +
  geom_boxplot() +
  labs(
    x = "Community cluster",
    y = "Mean warm-season CVD deaths per 100,000 (2011–2022)"
  ) +
  theme_classic()

p_box

ggsave(
  file.path(out_fig, "cluster_boxplot_CVD_rate.png"),
  p_box, width = 11, height = 8, dpi = 600
)


###### PCA + k-means for ED visits

df_ed_raw <- try_read_csv(paths$ed)
stopifnot(!is.null(df_ed_raw))

date_col <- c("ENC_ADMIT_DATE","ED_DATE","DATE")[c("ENC_ADMIT_DATE","ED_DATE","DATE") %in% names(df_ed_raw)][1]
stopifnot(!is.null(date_col))

df_ed <- df_ed_raw %>%
  mutate(DATE = as.Date(.data[[date_col]]))

df_ed <- ensure_community(df_ed)
if (!"community" %in% names(df_ed)) stop("Could not derive community for ED data")

present_ed <- intersect(expected_flags, names(df_ed))
if (!length(present_ed)) stop("No CVD/CHD/MI/Stroke flags found in ED file.")

df_ed <- df_ed %>%
  filter(month(DATE) %in% 5:9) %>%
  filter(year(DATE) >= cfg$start_year,
         year(DATE) <= cfg$end_year)

ed_daily <- df_ed %>%
  group_by(community, DATE) %>%
  summarise(
    across(all_of(present_ed), ~sum(.x, na.rm=TRUE), .names="{.col}_count"),
    all_visits = n(),
    .groups="drop"
  ) %>%
  mutate(
    year  = year(DATE),
    month = month(DATE),
    day   = day(DATE)
  )

ed_yearly <- ed_daily %>%
  group_by(community, year) %>%
  summarise(
    across(ends_with("_count"), sum, na.rm=TRUE),
    all_visits = sum(all_visits, na.rm=TRUE),
    .groups="drop"
  )

panel_ed_yearly <- acs_summary %>%
  inner_join(ed_yearly, by=c("community","year")) %>%
  inner_join(temps_yearly, by=c("community","year")) %>%
  mutate(all_deaths = all_visits) %>%
  left_join(ndvi_year, by=c("community","year")) %>%   # <-- DID YOU ADD THIS?
  mk_rates()

feat_ed <- panel_ed_yearly %>%
  filter(year >= 2011, year <= 2022) %>%
  group_by(community) %>%
  summarise(
    # ED outcomes
    CVD_ED = mean(CVD_per_100k, na.rm=TRUE),
    CHD_ED = mean(CHD_per_100k, na.rm=TRUE),
    MI_ED  = mean(MI_per_100k,  na.rm=TRUE),
    Stroke_ED = mean(Stroke_per_100k, na.rm=TRUE),
    
    # Climate
    mean_tmax     = mean(mean_tmax, na.rm=TRUE),
    mean_humidity = mean(mean_humidity, na.rm=TRUE),
    
    # ACS
    median_age     = mean(median_age, na.rm=TRUE),
    median_income  = mean(median_income, na.rm=TRUE),
    mean_college   = mean(mean_college, na.rm=TRUE),
    mean_hs        = mean(mean_hs, na.rm=TRUE),
    mean_white     = mean(mean_white, na.rm=TRUE),
    mean_black     = mean(mean_black, na.rm=TRUE),
    mean_hispanic  = mean(mean_hispanic, na.rm=TRUE),
    mean_asian     = mean(mean_asian, na.rm=TRUE),
    mean_male      = mean(mean_male, na.rm=TRUE),
    mean_female    = mean(mean_female, na.rm=TRUE),
    mean_employed  = mean(mean_employed, na.rm=TRUE),
    mean_unemployed= mean(mean_unemployed, na.rm=TRUE),
    
    # NDVI
    ndvi_mean      = mean(ndvi, na.rm=TRUE),
    
    .groups = "drop"
  )

# Matrix for PCA
mat <- feat_ed %>%
  column_to_rownames("community") %>%
  as.matrix()

# drop zero-variance columns
sd_vals <- apply(mat, 2, sd, na.rm = TRUE)
mat2 <- mat[, sd_vals > 0]

# scale
mat_scaled <- scale(mat2)

# PCA
pca_ed <- prcomp(mat_scaled)

# k-means
set.seed(42)
km_ed <- kmeans(pca_ed$x[,1:3], centers=3, nstart=100)

ed_clusters <- tibble(
  community = rownames(pca_ed$x),
  cluster = factor(km_ed$cluster)
)


# -----------------------------
# 1) Scree plot: variance explained (ED PCA)
# -----------------------------
imp_ed <- summary(pca_ed)$importance
var_exp_ed <- imp_ed[2, ]   # proportion of variance
cum_exp_ed <- imp_ed[3, ]   # cumulative

scree_ed <- tibble(
  PC   = factor(seq_along(var_exp_ed)),
  var  = var_exp_ed,
  cum  = cum_exp_ed
)

p_scree_ed <- ggplot(scree_ed, aes(x = as.numeric(PC), y = var)) +
  geom_col(width = 0.7) +
  geom_point() +
  geom_line(aes(y = cum), linetype = 2) +
  scale_x_continuous(
    breaks = seq_len(min(10, nrow(scree_ed))),
    labels = paste0("PC", seq_len(min(10, nrow(scree_ed))))
  ) +
  labs(
    x = NULL,
    y = "Proportion of variance explained",
    title = "ED PCA: Variance explained by principal components",
    subtitle = "Bars = individual PCs; dashed line = cumulative variance"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(out_fig, "ED_PCA_scree.png"),
  p_scree_ed, width = 7, height = 4.5, dpi = 600
)


# -----------------------------
# 2) Loadings table (PC1–PC3) for supplement
# -----------------------------
load_ed <- as.data.frame(pca_ed$rotation)
load_ed$variable <- rownames(load_ed)

loadings_ed_tbl <- load_ed %>%
  select(variable, PC1, PC2, PC3) %>%
  arrange(desc(abs(PC1)))   # or sort by any PC you care about

readr::write_csv(
  loadings_ed_tbl,
  file.path(out_tab, "ED_PCA_loadings_PC1_3.csv")
)


# -----------------------------
# 3) K-means clustering on ED PCA scores
#    (use first few PCs; adjust k as needed)
# -----------------------------
pc_scores_ed <- as.data.frame(pca_ed$x)
pc_scores_ed$community <- feat_ed$community

# how many PCs to use? often 2–4 is fine
pcs_use <- 1:3   # same as mortality, for consistency

set.seed(123)
k_ed <- 3        # change if elbow/silhouette suggests another k
km_ed <- kmeans(pc_scores_ed[, pcs_use], centers = k_ed, nstart = 100)

ed_clusters <- tibble(
  community = pc_scores_ed$community,
  cluster   = factor(km_ed$cluster)
)

readr::write_csv(
  ed_clusters,
  file.path(out_tab, "ED_clusters_assignments.csv")
)


# -----------------------------
# 4) PCA score plot (PC1 vs PC2) colored by ED cluster
# -----------------------------
scores_plot_ed <- pc_scores_ed %>%
  select(community, PC1, PC2) %>%
  left_join(ed_clusters, by = "community")

p_pca_ed <- ggplot(scores_plot_ed,
                   aes(x = PC1, y = PC2, color = cluster, label = community)) +
  geom_point(size = 2.4, alpha = 0.9) +
  # uncomment if you want labels:
  # ggrepel::geom_text_repel(size = 2.3, show.legend = FALSE) +
  labs(
    x = "PC1",
    y = "PC2",
    color = "Cluster",
    title = "Community clustering by ED-related heat burden and covariates",
    subtitle = "PCA scores (PC1 vs PC2), colored by k-means ED cluster"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_fig, "ED_PCA_PC1_PC2_clusters.png"),
  p_pca_ed, width = 7.2, height = 5.0, dpi = 600
)


# -----------------------------
# 5) ED cluster map (by community)
# -----------------------------
ed_map <- cas %>%
  left_join(ed_clusters, by = "community") %>%
  sf::st_as_sf()

tmap_mode("plot")

p_map_ed <- tm_shape(ed_map) +
  tm_polygons(
    "cluster",
    palette = "Set2",
    title = "ED cluster"
  ) +
  tm_layout(
    frame = FALSE,
    legend.outside = TRUE,
    legend.outside.position = "right",
    title = "Clusters of heat-associated ED burden and community characteristics"
  )

tmap::tmap_save(
  p_map_ed,
  filename = file.path(out_fig, "ED_clusters_map.png"),
  width = 1400, height = 1000, dpi = 200
)


# -----------------------------
# 6) Cluster summary table (for text / supplement)
# -----------------------------
cluster_summary_ed <- feat_ed %>%
  left_join(ed_clusters, by = "community") %>%
  group_by(cluster) %>%
  summarise(
    n_communities = dplyr::n(),
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

readr::write_csv(
  cluster_summary_ed,
  file.path(out_tab, "ED_cluster_summary.csv")
)

# Merge ED CVD rate with clusters
ed_plot_dat <- feat_ed %>%
  left_join(ed_clusters, by = "community")

# Quick check
# head(ed_plot_dat)

p_box_ed <- ggplot(ed_plot_dat, aes(x = cluster, y = CVD_ED)) +
  geom_boxplot(outlier.alpha = 0.4, width = 0.6) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 1.8) +
  labs(
    x = "Community vulnerability cluster",
    y = "CVD ED visit rate (per 100,000; 2011–2022 warm seasons)",
    title = "Distribution of heat-season CVD ED visit rates by community cluster"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_box_ed

ggsave(
  file.path(out_fig, "ED_CVD_rate_by_cluster_boxplot.png"),
  p_box_ed, width = 6.5, height = 4.8, dpi = 600
)

ggsave(
  file.path(out_fig, "ED_CVD_rate_by_cluster_boxplot.pdf"),
  p_box_ed, width = 6.5, height = 4.8, device = cairo_pdf
)




imp_ed <- summary(pca_ed)$importance
var_exp_ed <- imp_ed[2, ]   # proportion of variance
cum_exp_ed <- imp_ed[3, ]   # cumulative proportion

scree_ed <- tibble(
  PC  = seq_along(var_exp_ed),
  var = var_exp_ed,
  cum = cum_exp_ed
)

p_scree_ed <- ggplot(scree_ed, aes(x = PC)) +
  geom_col(aes(y = var), width = 0.7) +
  geom_point(aes(y = cum)) +
  geom_line(aes(y = cum), linetype = 2) +
  scale_x_continuous(
    breaks = scree_ed$PC,
    labels = paste0("PC", scree_ed$PC)
  ) +
  labs(
    x = NULL,
    y = "Proportion of variance explained",
    title = "Principal component variance explained (ED analysis)",
    subtitle = "Bars: individual components; dashed line: cumulative variance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    plot.subtitle= element_text(hjust = 0.5)
  )

p_scree_ed

ggsave(
  file.path(out_fig, "ED_PCA_scree.png"),
  p_scree_ed, width = 7.0, height = 4.8, dpi = 600
)

ggsave(
  file.path(out_fig, "ED_PCA_scree.pdf"),
  p_scree_ed, width = 7.0, height = 4.8, device = cairo_pdf
)

library(ggrepel)

scores_ed <- as.data.frame(pca_ed$x)
scores_ed$community <- feat_ed$community

scores_ed <- scores_ed %>%
  left_join(ed_clusters, by = "community")

p_pca_ed <- ggplot(scores_ed,
                   aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.4, alpha = 0.9) +
  # uncomment to label communities:
  ggrepel::geom_text_repel(aes(label = community), size = 2.3, show.legend = FALSE) +
  labs(
    x = "PC1",
    y = "PC2",
    color = "Cluster",
    title = "Community vulnerability clusters in PCA space",
    subtitle = "ED analysis: PC1 vs PC2 scores colored by k-means cluster"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

p_pca_ed

ggsave(
  file.path(out_fig, "ED_PCA_PC1_PC2_clusters.png"),
  p_pca_ed, width = 7.2, height = 5.0, dpi = 600
)

ggsave(
  file.path(out_fig, "ED_PCA_PC1_PC2_clusters.pdf"),
  p_pca_ed, width = 7.2, height = 5.0, device = cairo_pdf
)

load_ed <- as.data.frame(pca_ed$rotation)
load_ed$variable <- rownames(load_ed)

loadings_ed_tbl <- load_ed %>%
  select(variable, PC1, PC2, PC3) %>%
  arrange(desc(abs(PC1)))   # or sort by whichever PC you talk about most

readr::write_csv(
  loadings_ed_tbl,
  file.path(out_tab, "ED_PCA_loadings_PC1_3.csv")
)

cluster_summary_ed <- feat_ed %>%
  left_join(ed_clusters, by = "community") %>%
  group_by(cluster) %>%
  summarise(
    n_communities = dplyr::n(),
    CVD_ED_mean   = mean(CVD_ED, na.rm = TRUE),
    CVD_ED_sd     = sd(CVD_ED, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  cluster_summary_ed,
  file.path(out_tab, "ED_CVD_rate_by_cluster_summary.csv")
)

# Merge mortality + ED PCA scores and clusters

mort_scores <- feat %>%
  dplyr::select(community) %>%
  dplyr::mutate(
    PC1_mort = pca$x[, 1],
    PC2_mort = pca$x[, 2]
  ) %>%
  dplyr::left_join(
    clusters %>% dplyr::select(community, cluster_mort = cluster),
    by = "community"
  )

ed_scores <- feat_ed %>%
  dplyr::select(community) %>%
  dplyr::mutate(
    PC1_ed = pca_ed$x[, 1],
    PC2_ed = pca_ed$x[, 2]
  ) %>%
  dplyr::left_join(
    ed_clusters %>% dplyr::select(community, cluster_ed = cluster),
    by = "community"
  )

compare_df <- mort_scores %>%
  dplyr::inner_join(ed_scores, by = "community")


# Make sure cluster variables are factors with nice labels
compare_df <- compare_df %>%
  mutate(
    cluster_mort = factor(cluster_mort),
    cluster_ed   = factor(cluster_ed)
  )

p_vuln_scatter <- ggplot(
  compare_df,
  aes(x = PC1_mort, y = PC1_ed,
      color = cluster_mort, shape = cluster_ed, label = community)
) +
  # reference lines
  geom_hline(yintercept = 0, linetype = 3, color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = 3, color = "grey70", linewidth = 0.4) +
  # 1:1 line: equal mortality & ED vulnerability
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "black", linewidth = 0.5) +
  
  # points + optional labels
  geom_point(size = 2.8, alpha = 0.9) +
  ggrepel::geom_text_repel(
    size = 2.4,
    max.overlaps = 50,
    segment.color = "grey70",
    show.legend = FALSE
  ) +
  
  scale_color_brewer(palette = "Dark2", name = "Mortality cluster") +
  scale_shape_manual(
    values = c(16, 17, 15, 3)[seq_along(levels(compare_df$cluster_ed))],
    name   = "ED visit cluster"
  ) +
  
  labs(
    x = "CVD mortality vulnerability (PC1 score)",
    y = "CVD ED visit vulnerability (PC1 score)",
    title = "Comparing community vulnerability to heat-related CVD mortality vs ED visits",
    subtitle = "Points are Chicago community areas; colors = mortality clusters, shapes = ED clusters\nDashed line shows equal mortality and ED vulnerability"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    legend.box      = "vertical"
  )

p_vuln_scatter

ggsave(
  file.path(out_fig, "fig_mort_vs_ed_vulnerability_scatter.png"),
  p_vuln_scatter, width = 7.5, height = 5.6, dpi = 600
)

ggsave(
  file.path(out_fig, "fig_mort_vs_ed_vulnerability_scatter.pdf"),
  p_vuln_scatter, width = 7.5, height = 5.6, device = cairo_pdf
)


# -------------------------------
# 2) Cross-classification heatmap: mortality × ED clusters
# -------------------------------

cross_tab <- compare_df %>%
  count(cluster_mort, cluster_ed, name = "n_communities")

p_cross_heat <- ggplot(
  cross_tab,
  aes(x = cluster_ed, y = cluster_mort, fill = n_communities)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = n_communities), size = 4.0) +
  scale_fill_gradient(
    low  = "white",
    high = "steelblue",
    name = "Number of\ncommunities",
    na.value = "grey90"
  ) +
  labs(
    x = "ED visit vulnerability cluster",
    y = "Mortality vulnerability cluster",
    title = "Cross-classification of community vulnerability",
    subtitle = "Counts of Chicago community areas by joint mortality and ED vulnerability clusters"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

p_cross_heat

ggsave(
  file.path(out_fig, "fig_mort_vs_ed_cluster_crosstab_heatmap.png"),
  p_cross_heat, width = 6.4, height = 4.8, dpi = 600
)

ggsave(
  file.path(out_fig, "fig_mort_vs_ed_cluster_crosstab_heatmap.pdf"),
  p_cross_heat, width = 6.4, height = 4.8, device = cairo_pdf
)

# also save the cross-tab as a table for the supplement
readr::write_csv(
  cross_tab,
  file.path(out_tab, "mort_vs_ed_cluster_crosstab.csv")
)


# --- Combined dataset: mortality + ED + cluster ---
mort_rates <- panel_yearly_simple %>%
  filter(year >= 2011, year <= 2022) %>%
  group_by(community) %>%
  summarise(
    CVD_mort = mean(CVD_per_100k, na.rm = TRUE),
    .groups = "drop"
  )

ed_rates <- panel_ed_yearly %>%
  filter(year >= 2011, year <= 2022) %>%
  group_by(community) %>%
  summarise(
    CVD_ED = mean(CVD_per_100k, na.rm = TRUE),
    .groups = "drop"
  )

comm_rates <- mort_rates %>%
  inner_join(ed_rates, by = "community") %>%
  inner_join(clusters, by = "community")

p_scatter <- ggplot(
  comm_rates,
  aes(x = CVD_mort, y = CVD_ED,
      color = cluster, label = community)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 60) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 2, color = "black") +
  labs(
    x = "CVD Mortality Rate (per 100k)",
    y = "CVD ED Visit Rate (per 100k)",
    title = "Community-level CVD Mortality vs ED Visit Rates",
    subtitle = "Color indicates PCA-derived vulnerability cluster"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic(base_size = 18)

p_scatter

rates_long <- comm_rates %>%
  tidyr::pivot_longer(
    cols = c(CVD_mort, CVD_ED),
    names_to = "outcome",
    values_to = "rate"
  ) %>%
  mutate(
    outcome = recode(outcome,
                     CVD_mort = "CVD Mortality",
                     CVD_ED   = "CVD ED Visits")
  )

p_box <- ggplot(
  rates_long,
  aes(x = cluster, y = rate, fill = outcome)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.55,
    width = 0.55,
    position = position_dodge(width = 0.65)
  ) +
  geom_jitter(
    aes(color = outcome),
    size = 2.0,
    alpha = 0.85,
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.65
    )
  ) +
  scale_fill_brewer(palette = "Set1", name = NULL) +     # ← no legend title
  scale_color_brewer(palette = "Set1", guide = "none") +
  coord_cartesian(ylim = c(0, 1500)) +                   # ← y-axis cap
  labs(
    x = "Vulnerability cluster",
    y = "Rate per 100,000 population",
    title = "CVD Mortality and ED Visit Rates by Vulnerability Cluster"
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 16),               # ← larger legend text
    legend.key.size = unit(1.9, "lines"),                # ← larger legend key
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_box

# ---- SAVE SCATTER ----
ggsave(
  filename = file.path(out_fig, "scatter_CVD_mort_vs_ED_by_cluster.png"),
  plot     = p_scatter,
  width    = 11.2, height = 9.2, dpi = 600
)

ggsave(
  filename = file.path(out_fig, "scatter_CVD_mort_vs_ED_by_cluster.pdf"),
  plot     = p_scatter,
  width    = 11.2, height = 9.2, device = cairo_pdf
)

# ---- SAVE BOXPLOT ----
ggsave(
  filename = file.path(out_fig, "box_CVD_mort_ED_by_cluster.png"),
  plot     = p_box,
  width    = 11.2, height = 9.2, dpi = 600
)

ggsave(
  filename = file.path(out_fig, "box_CVD_mort_ED_by_cluster.pdf"),
  plot     = p_box,
  width    = 11.2, height = 9.2, device = cairo_pdf
)


# --- join cluster assignment to shapefile ---
cas_cluster <- cas %>%
  left_join(comm_rates %>% select(community, cluster), by = "community") %>%
  mutate(cluster = factor(cluster))

# --- map ---
p_map <- ggplot(cas_cluster) +
  geom_sf(aes(fill = cluster), color = "white", size = 0.35) +
  scale_fill_brewer(palette = "Dark2", name = "Cluster") +
  labs(
    title = "Community Vulnerability Clusters",
    subtitle = "Derived from PCA of heat exposure, demographics, NDVI, and air pollution"
    # caption removed
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 14),         # ← bigger legend text
    legend.title = element_text(size = 13, face = "bold"),
    legend.key.size = unit(3.0, "lines"),          # ← bigger keys
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

p_map

ggsave(
  filename = file.path(out_fig, "map_clusters_CVD_vulnerability.png"),
  plot     = p_map,
  width    = 11, height = 13, dpi = 600
)

ggsave(
  filename = file.path(out_fig, "map_clusters_CVD_vulnerability.pdf"),
  plot     = p_map,
  width    = 11, height = 13, device = cairo_pdf
)

# --- Scree data ---
eig_vals <- pca$sdev^2
var_explained <- eig_vals / sum(eig_vals)

scree_df <- tibble(
  PC        = factor(seq_along(var_explained)),
  var_prop  = var_explained,
  cum_prop  = cumsum(var_explained)
)

# --- Scree plot: bars + cumulative line ---
p_scree <- ggplot(scree_df, aes(x = PC, y = var_prop)) +
  geom_col(fill = "grey70", color = "grey40") +
  geom_point(aes(y = cum_prop), size = 2.2, color = "steelblue4", group = 1) +
  geom_line(aes(y = cum_prop), linewidth = 0.7, color = "steelblue4", group = 1) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Principal component",
    y = "Variance explained",
    title = "Scree plot of principal components",
    subtitle = "Bars: proportion of variance explained; line: cumulative"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

p_scree

ggsave(
  file.path(out_fig, "pca_scree_plot.png"),
  p_scree, width = 6.5, height = 4.5, dpi = 600
)

ggsave(
  file.path(out_fig, "pca_scree_plot.pdf"),
  p_scree, width = 6.5, height = 4.5, device = cairo_pdf
)


# --- Scores with community + cluster ---
scores_df <- as.data.frame(pca$x[, 1:2]) %>%
  mutate(
    community = feat$community  # assumes same row order as pca input
  ) %>%
  left_join(
    clusters,                        # has community, cluster
    by = "community"
  ) %>%
  mutate(cluster = factor(cluster))

# rename PC columns for clarity
scores_df <- scores_df %>%
  rename(
    PC1 = PC1,
    PC2 = PC2
  )

# --- PC1 vs PC2 plot ---
p_pc12 <- ggplot(scores_df,
                 aes(x = PC1, y = PC2,
                     color = cluster, label = community)) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey80", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = 3, color = "grey80", linewidth = 0.4) +
  geom_point(size = 2.8, alpha = 0.9) +
  ggrepel::geom_text_repel(
    size = 2.5,
    max.overlaps = 60,
    segment.color = "grey75",
    show.legend = FALSE
  ) +
  scale_color_brewer(palette = "Dark2", name = "Cluster") +  # same palette as map
  labs(
    x = "PC1",
    y = "PC2",
    title = "PCA of community vulnerability features",
    subtitle = "PC1 vs PC2 scores for Chicago community areas, colored by cluster"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

p_pc12

ggsave(
  file.path(out_fig, "pca_PC1_PC2_clusters.png"),
  p_pc12, width = 7.0, height = 5.3, dpi = 600
)

ggsave(
  file.path(out_fig, "pca_PC1_PC2_clusters.pdf"),
  p_pc12, width = 7.0, height = 5.3, device = cairo_pdf
)


# Loadings as a data frame
loadings <- as.data.frame(pca$rotation[, 1:2])
loadings$variable <- rownames(loadings)

# Top contributors to PC1
top_pc1 <- loadings %>%
  arrange(desc(abs(PC1))) %>%
  select(variable, PC1) %>%
  slice_head(n = 10)

# Top contributors to PC2
top_pc2 <- loadings %>%
  arrange(desc(abs(PC2))) %>%
  select(variable, PC2) %>%
  slice_head(n = 10)

top_pc1
top_pc2

pc1_label <- "PC1: socioeconomic disadvantage & heat exposure\n(↑ unemployment, ↑ heat, ↑ Black; ↓ income, ↓ college)"
pc2_label <- "PC2: racial/ethnic composition & humidity contrast\n(↓ Hispanic vs ↑ Black/humid communities)"

p_pc12 <- ggplot(scores_df,
                 aes(x = PC1, y = PC2,
                     color = cluster, label = community)) +
  geom_hline(yintercept = 0, linetype = 3, color = "grey80", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = 3, color = "grey80", linewidth = 0.4) +
  geom_point(size = 2.8, alpha = 0.9) +
  ggrepel::geom_text_repel(
    size = 2.5,
    max.overlaps = 60,
    segment.color = "grey75",
    show.legend = FALSE
  ) +
  scale_color_brewer(palette = "Dark2", name = "Cluster") +
  labs(
    x = pc1_label,
    y = pc2_label,
    title = "PCA of Community Vulnerability Features",
    subtitle = "PC1 vs PC2 scores for Chicago community areas"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

p_pc12

ggsave(
  file.path(out_fig, "pca_PC1_PC2_clusters.png"),
  p_pc12, width = 12.0, height = 8.3, dpi = 600
)


loadings_full <- as.data.frame(pca$rotation[, 1:2])
loadings_full$variable <- rownames(loadings_full)

loadings_full <- loadings_full %>%
  relocate(variable, PC1, PC2) %>%
  arrange(desc(abs(PC1)))  # or abs(PC2) if you prefer

loadings_full


### maps of NDVI and NO2 by community area

library(terra)       # for reading NetCDF and spatial raster ops
library(sf)          # for vector geometries
library(dplyr)
library(ggplot2)
library(exactextractr)  # fast zonal extraction



ndvi_dir <- "/Users/saborpete/Desktop/Peter/PhD/Dissertation/Analysis/Code and Data/Heat_CVD_Chicago/data/ndvi/NCEI"

years <- 2011:2022

ndvi_yearly <- list()

for (yr in years) {
  
  nc_path <- list.files(file.path(ndvi_dir, as.character(yr)),
                        pattern = "\\.nc$", full.names = TRUE)
  
  if (length(nc_path) == 0) next
  
  message("Reading NDVI year: ", yr)
  
  r <- rast(nc_path)
  
  # crop raster to Chicago extent
  r_chi <- crop(r, vect(cas))
  r_chi <- mask(r_chi, vect(cas))
  
  # compute community-area mean NDVI
  ndvi_vals <- exact_extract(r_chi, cas, "mean")
  
  ndvi_yearly[[as.character(yr)]] <- 
    tibble(
      community = cas$community,
      year = yr,
      ndvi = ndvi_vals
    )
}

ndvi_df <- bind_rows(ndvi_yearly)

# multi-year mean NDVI for mapping
ndvi_mean <- ndvi_df %>%
  group_by(community) %>%
  summarize(ndvi_mean = mean(ndvi$mean.QA, na.rm = TRUE))


cas_ndvi <- cas %>%
  left_join(ndvi_mean, by = "community")

p_ndvi <- ggplot(cas_ndvi) +
  geom_sf(aes(fill = ndvi_mean), color = NA) +
  scale_fill_viridis_c(option = "viridis", name = "Mean NDVI") +
  labs(
    title = "Mean 2011–2022 NDVI",
    subtitle = "Averaged across MODIS/AVHRR NOAA CDR NDVI dataset"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_ndvi

ggsave(
  filename = file.path(out_fig, "map_ndvi_2011_2022.png"),
  plot = p_ndvi,
  width = 7.5, height = 9, dpi = 600
)

no2_path <- "/Users/saborpete/Desktop/Peter/PhD/Dissertation/Analysis/Code and Data/Heat_CVD_Chicago/data/2019_final_1km.nc"

# read raster, crop, mask
r_no2 <- rast(no2_path)

r_no2_chi <- crop(r_no2, vect(cas))
r_no2_chi <- mask(r_no2_chi, vect(cas))

# extract community means
no2_vals <- exact_extract(r_no2_chi, cas, "mean")

no2_df <- tibble(
  community = cas$community,
  no2 = no2_vals
)

cas_no2 <- cas %>%
  left_join(no2_df, by = "community")

p_no2 <- ggplot(cas_no2) +
  geom_sf(aes(fill = no2), color = NA) +
  scale_fill_viridis_c(option = "inferno", name = "NO2 (µg/m³)") +
  labs(
    title = "Mean NO2 Concentrations",
    subtitle = "2019 Global 1-km Dataset"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_no2


library(terra)
library(sf)
library(dplyr)
library(stringr)
library(exactextractr)
library(ggplot2)

# directories you gave
pm25_dir <- "/Users/saborpete/Desktop/Peter/Postdoc/Data/Annual-selected"
no2_dir  <- "/Users/saborpete/Desktop/Peter/Postdoc/Data/Anenberg et al 2022 Global NO2-selected"


process_pollutant_nc <- function(folder, value_name = "value") {
  files <- list.files(folder, pattern = "\\.nc$", full.names = TRUE)
  if (!length(files)) stop("No .nc files found in: ", folder)
  
  all_years <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    f <- files[i]
    base <- basename(f)
    
    # extract 4-digit year from filename
    year <- str_extract(base, "20[0-9]{2}")
    if (is.na(year)) {
      warning("Could not detect year in filename: ", base, "; skipping.")
      next
    }
    year <- as.integer(year)
    
    message("Processing ", base, " (year = ", year, ")")
    
    r <- rast(f)
    # if multiple layers (e.g., monthly slices), average them first
    if (nlyr(r) > 1) {
      r <- app(r, mean, na.rm = TRUE)
    }
    
    # crop + mask to Chicago
    r_chi <- crop(r, vect(cas))
    r_chi <- mask(r_chi, vect(cas))
    
    # mean value per community area
    vals <- exact_extract(r_chi, cas, "mean")
    
    all_years[[i]] <- tibble(
      community = cas$community,
      year      = year,
      !!value_name := vals
    )
  }
  
  bind_rows(all_years) %>%
    arrange(community, year)
}

pm25_panel <- process_pollutant_nc(pm25_dir, value_name = "pm25")
no2_panel  <- process_pollutant_nc(no2_dir,  value_name = "no2")

# optional: save for your modeling pipeline
write.csv(pm25_panel,
          "/Users/saborpete/Desktop/Peter/PhD/Dissertation/Analysis/Code and Data/Heat_CVD_Chicago/output/tables/pm25_panel_community_year.csv",
          row.names = FALSE)

write.csv(no2_panel,
          "/Users/saborpete/Desktop/Peter/PhD/Dissertation/Analysis/Code and Data/Heat_CVD_Chicago/output/tables/no2_panel_community_year.csv",
          row.names = FALSE)



pm25_mean <- pm25_panel %>%
  group_by(community) %>%
  summarise(pm25_mean = mean(pm25, na.rm = TRUE), .groups = "drop")

no2_mean <- no2_panel %>%
  group_by(community) %>%
  summarise(no2_mean = mean(no2, na.rm = TRUE), .groups = "drop")


cas_pm25 <- cas %>% left_join(pm25_mean, by = "community")
cas_no2  <- cas %>% left_join(no2_mean,  by = "community")


p_pm25 <- ggplot(cas_pm25) +
  geom_sf(aes(fill = pm25_mean), color = "white", size = 0.2) +
  scale_fill_distiller(
    palette = "PuBuGn",
    direction = 1,
    name = bquote("Mean PM"[2.5] * " (µg/m³)")
  ) +
  labs(
    title = bquote("Mean 2011–2022 PM"[2.5]),
    subtitle = "Annual average fine particulate matter by Chicago community area"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

p_pm25

p_no2 <- ggplot(cas_no2) +
  geom_sf(aes(fill = no2_mean), color = "white", size = 0.2) +
  scale_fill_distiller(
    palette = "OrRd",
    direction = 1,
    name = bquote("Mean NO"[2] * " (ppb)")
  ) +
  labs(
    title = bquote("Mean 2011–2022 NO"[2]),
    subtitle = "Annual average nitrogen dioxide by Chicago community area"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )


p_no2

# Ensure output folder exists
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

# ---- Save PM2.5 map ----
ggsave(
  filename = file.path(out_fig, "map_pm25_2011_2022.png"),
  plot = p_pm25,
  width = 7.5, height = 9, dpi = 600
)

ggsave(
  filename = file.path(out_fig, "map_pm25_2011_2022.pdf"),
  plot = p_pm25,
  width = 7.5, height = 9, device = cairo_pdf
)

# ---- Save NO2 map ----
ggsave(
  filename = file.path(out_fig, "map_no2_2011_2022.png"),
  plot = p_no2,
  width = 7.5, height = 9, dpi = 600
)

ggsave(
  filename = file.path(out_fig, "map_no2_2011_2022.pdf"),
  plot = p_no2,
  width = 7.5, height = 9, device = cairo_pdf
)


full_panel <- temps_yearly %>%
  left_join(ndvi_year, by = c("community","year")) %>%
  left_join(no2_panel,  by = c("community","year")) %>%
  left_join(pm25_panel, by = c("community","year")) %>%
  left_join(acs_summary, by = c("community","year"))


full_panel <- full_panel %>%
  filter(year >= 2011, year <= 2022) %>%
  group_by(community) %>%
  summarise(
    mean_tmax = mean(mean_tmax, na.rm = TRUE),
    ndvi      = mean(ndvi, na.rm = TRUE),
    no2       = mean(no2, na.rm = TRUE),
    pm25      = mean(pm25, na.rm = TRUE),
    
    # ACS
    median_age     = mean(median_age, na.rm = TRUE),
    median_income  = mean(median_income, na.rm = TRUE),
    mean_college   = mean(mean_college, na.rm = TRUE),
    mean_hs        = mean(mean_hs, na.rm = TRUE),
    mean_white     = mean(mean_white, na.rm = TRUE),
    mean_black     = mean(mean_black, na.rm = TRUE),
    mean_hispanic  = mean(mean_hispanic, na.rm = TRUE),
    mean_asian     = mean(mean_asian, na.rm = TRUE),
    mean_employed  = mean(mean_employed, na.rm = TRUE),
    mean_unemployed= mean(mean_unemployed, na.rm = TRUE),
    mean_male      = mean(mean_male, na.rm = TRUE),
    mean_female    = mean(mean_female, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    tmax_tertile = ntile(mean_tmax, 3)
  )

summary_median_iqr <- function(x) {
  med <- median(x, na.rm = TRUE)
  q1  <- quantile(x, 0.25, na.rm = TRUE)
  q3  <- quantile(x, 0.75, na.rm = TRUE)
  sprintf("%.1f (%.1f–%.1f)", med, q1, q3)
}

table2_env <- full_panel %>%
  group_by(tmax_tertile) %>%
  summarise(
    NDVI         = summary_median_iqr(ndvi),
    NO2          = summary_median_iqr(no2),
    PM25         = summary_median_iqr(pm25),
    
    .groups = "drop"
  ) %>%
  mutate(
    tmax_tertile = c("Low Tmax Tertile", "Medium Tmax Tertile", "High Tmax Tertile")
  )


# ---- paths ----
root <- "/Users/saborpete/Desktop/Peter/PhD/Dissertation/Analysis/Code and Data/Heat_CVD_Chicago"
tab_dir <- file.path(root, "output", "tables")

# 1) Load the peak-threshold / excess-rate tables you already created
results_ED   <- read_csv(file.path(tab_dir, "results_ED.csv"),   show_col_types = FALSE)
results_MORT <- read_csv(file.path(tab_dir, "results_MORT.csv"), show_col_types = FALSE)

# 2) Load the GAM edf/p summaries
gam_ED   <- read_csv(file.path(tab_dir, "gam_summary_ED_edf_p_r2_dev.csv"),
                     show_col_types = FALSE)
gam_MORT <- read_csv(file.path(tab_dir, "gam_summary_MORT_edf_p_r2_dev.csv"),
                     show_col_types = FALSE)

# 3) Make temporal labels line up (annual ↔ Warm Months, etc.)
recode_temporal <- function(x) {
  recode(x,
         "yearly"  = "annual",
         "monthly" = "monthly",
         "daily"   = "daily",
         "dailyMA" = "dailyMA")
}

gam_ED2 <- gam_ED %>%
  mutate(temporal = recode_temporal(scale)) %>%
  transmute(
    subtype  = outcome,
    temporal,
    p_value  = p
  )

gam_MORT2 <- gam_MORT %>%
  mutate(temporal = recode_temporal(scale)) %>%
  transmute(
    subtype  = outcome,
    temporal,
    p_value  = p
  )

# 4) Join p-values onto the results tables
results_ED_p <- results_ED %>%
  left_join(gam_ED2, by = c("subtype", "temporal"))

results_MORT_p <- results_MORT %>%
  left_join(gam_MORT2, by = c("subtype", "temporal"))

# 5) Nicely formatted p-value column for the manuscript
format_p <- function(x) {
  if (is.na(x)) return(NA_character_)
  if (x < 0.001) return("<0.001")
  if (x < 0.01)  return(sprintf("%.3f", x))
  sprintf("%.3f", x)
}

results_ED_p  <- results_ED_p  %>% mutate(p_display = vapply(p_value, format_p, character(1)))
results_MORT_p<- results_MORT_p%>% mutate(p_display = vapply(p_value, format_p, character(1)))

# 6) If you want a single Table 2-style object:
table2 <- bind_rows(
  results_ED_p  %>% mutate(outcome_type = "Emergency Department Visits"),
  results_MORT_p%>% mutate(outcome_type = "Mortality")
) %>%
  mutate(
    temporal_label = recode(temporal,
                            "annual"  = "Warm Months",
                            "daily"   = "Daily",
                            "dailyMA" = "Daily (lag 0–3)",
                            "monthly" = "Monthly"
    )
  ) %>%
  select(
    `Outcome Type`        = outcome_type,
    `CVD Subtype`         = subtype,
    `Temporal Scale`      = temporal_label,
    `Peak Tmax Threshold` = peak_threshold,
    `Peak Excess Rate`    = peak_excess_rate,
    `Peak SD`             = peak_sd,
    `p-value`             = p_display
  )

# 7) Write out for copy/paste into Word
write_csv(table2, file.path(tab_dir, "Table2_with_pvalues.csv"))


library(mgcv)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(gratia)

library(dplyr)
library(purrr)
library(stringr)
library(gratia)

extract_importance <- function(obj, model_name) {
  
  # unwrap if saved via save_gam_with_means()
  fit <- if (is.list(obj) && "fit" %in% names(obj)) obj$fit else obj
  
  sm <- summary(fit)
  
  ## ---- Smooth terms ----
  s_tab   <- as.data.frame(sm$s.table)
  s_terms <- rownames(sm$s.table)
  
  # Very safe spline direction function
  get_spline_direction <- function(term) {
    d <- tryCatch(
      gratia::derivatives(fit, select = term),
      error = function(e) NULL
    )
    if (is.null(d)) return(NA_character_)
    if (!("derivative" %in% names(d))) return(NA_character_)
    if (!is.numeric(d$derivative)) return(NA_character_)
    if (all(is.na(d$derivative))) return(NA_character_)
    
    mean_deriv <- mean(d$derivative, na.rm = TRUE)
    if (is.nan(mean_deriv)) return(NA_character_)
    if (mean_deriv > 0)  "positive"
    else if (mean_deriv < 0) "negative"
    else NA_character_
  }
  
  spline_dirs <- if (nrow(s_tab)) {
    purrr::map_chr(s_terms, get_spline_direction)
  } else character(0)
  
  smooth_df <- if (nrow(s_tab)) {
    tibble::tibble(
      model     = model_name,
      term      = s_terms,
      type      = "smooth",
      EDF       = s_tab$edf,
      F         = if ("F" %in% names(s_tab)) s_tab$F else NA_real_,
      p         = s_tab$`p-value`,
      direction = spline_dirs
    )
  } else tibble::tibble()
  
  ## ---- Parametric / linear terms ----
  p_tab   <- as.data.frame(sm$p.table)
  p_terms <- rownames(sm$p.table)
  
  # figure out which columns are the stat and p-value
  stat_col <- intersect(c("t value", "z value"), colnames(p_tab))
  stat_col <- if (length(stat_col)) stat_col[1] else NA_character_
  
  p_col <- grep("^Pr\\(", colnames(p_tab), value = TRUE)
  p_col <- if (length(p_col)) p_col[1] else NA_character_
  
  param_df <- if (nrow(p_tab)) {
    tibble::tibble(
      model     = model_name,
      term      = p_terms,
      type      = "linear",
      estimate  = p_tab$Estimate,
      se        = p_tab$Std..Error,
      stat      = if (!is.na(stat_col)) p_tab[[stat_col]] else NA_real_,
      p         = if (!is.na(p_col))    p_tab[[p_col]]    else NA_real_,
      direction = ifelse(p_tab$Estimate > 0, "positive", "negative")
    )
  } else tibble::tibble()
  
  dplyr::bind_rows(smooth_df, param_df)
}




model_dir   <- "/Users/saborpete/Desktop/Peter/PhD/Dissertation/Analysis/Code and Data/Heat_CVD_Chicago/output/models"
model_files <- list.files(model_dir, full.names = TRUE, pattern = "\\.rds$")

importance_tbl <- purrr::map_dfr(model_files, function(path) {
  model_name <- basename(path) %>% str_remove("\\.rds$")
  obj <- readRDS(path)
  extract_importance(obj, model_name)
})

# quick sanity check
importance_tbl %>%
  dplyr::filter(type == "linear") %>%
  dplyr::count(is_na_p = is.na(p))

# save for supplement
write_csv(
  importance_tbl,
  file.path(model_dir, "Supplement_Table_GAM_FeatureImportance.csv")
)

# Filter to only statistically significant terms
sig_tbl <- importance_tbl %>%
  filter(!is.na(p) & p < 0.05) %>%
  arrange(model, type, term)

# Write out for supplement
write_csv(
  sig_tbl,
  file.path(model_dir, "Supplement_Table_GAM_SignificantTerms.csv")
)

# Keep only daily models
sig_daily <- sig_tbl %>%
  filter(str_detect(model, "_daily\\.rds$") |      # strict pattern
           str_detect(model, "_daily$") |            # in case extension was stripped
           str_detect(model, "daily_")) %>%          # catch all forms
  arrange(model, type, term)

# Export
write_csv(
  sig_daily,
  file.path(model_dir, "Supplement_Table_GAM_SignificantTerms_Daily.csv")
)

sig_daily


library(DiagrammeR)

grViz("
digraph flowchart {
  graph [layout = dot, rankdir = TB]

  node [shape = box, fontsize = 12, width = 3]

  deaths_raw       [label = 'deaths_raw\\n1993–2022']
  yr               [label = 'Restrict to 2011–2022']
  ca               [label = 'Valid Chicago CA_NAME']
  warm             [label = 'Warm months\\n(May–September)']
  cvd_split        [label = 'CVD vs non-CVD']
  cvd_only         [label = 'CVD only']
  subtype          [label = 'CHD / MI / Stroke']
  merge            [label = 'Merge with climate, ACS,\nNDVI, PM2.5, NO2']

  deaths_raw -> yr -> ca -> warm -> cvd_split
  cvd_split -> cvd_only
  cvd_only -> subtype -> merge
}
")

make_flow_table <- function(df,
                            date_col,
                            community_col,
                            cvd_col   = "CVD",
                            chd_col   = "CHD",
                            mi_col    = "MI",
                            stroke_col= "Stroke",
                            start_year = 2011,
                            end_year   = 2022,
                            label_prefix = "Mortality") {
  
  df <- df %>%
    mutate(
      DATE = as.Date(.data[[date_col]]),
      year = year(DATE),
      month = month(DATE)
    )
  
  # 1) Raw rows
  n_raw <- nrow(df)
  
  # 2) Restrict to study years
  df_year <- df %>% filter(year >= start_year, year <= end_year)
  n_year <- nrow(df_year)
  
  # 3) Valid Chicago community
  df_comm <- df_year %>%
    filter(!is.na(.data[[community_col]]) & .data[[community_col]] != "")
  n_comm <- nrow(df_comm)
  
  # 4) Warm months (May–September)
  df_warm <- df_comm %>%
    filter(month %in% 5:9)
  n_warm <- nrow(df_warm)
  
  # 5) CVD vs non-CVD
  has_cvd <- df_warm[[cvd_col]] %in% c(1, TRUE)
  n_cvd   <- sum(has_cvd, na.rm = TRUE)
  n_non   <- n_warm - n_cvd
  
  # 6) Subtypes (among warm-month records, counted separately)
  n_chd    <- sum(df_warm[[chd_col]]    %in% c(1, TRUE), na.rm = TRUE)
  n_mi     <- sum(df_warm[[mi_col]]     %in% c(1, TRUE), na.rm = TRUE)
  n_stroke <- sum(df_warm[[stroke_col]] %in% c(1, TRUE), na.rm = TRUE)
  
  tibble::tibble(
    Dataset   = label_prefix,
    Step      = c(
      "Raw records",
      sprintf("Records %d–%d", start_year, end_year),
      "With valid Chicago community area",
      "Warm months only (May–September)",
      "CVD events (any CVD flag = 1)",
      "Non-CVD events",
      "CHD events",
      "MI events",
      "Stroke events"
    ),
    N         = c(
      n_raw,
      n_year,
      n_comm,
      n_warm,
      n_cvd,
      n_non,
      n_chd,
      n_mi,
      n_stroke
    )
  )
}

flow_mort <- make_flow_table(
  df            = deaths_raw,
  date_col      = "DOD",      # change if your date column is named differently
  community_col = "CA_NAME",  # change if needed
  label_prefix  = "Mortality"
)

flow_mort


flow_ed <- make_flow_table(
  df            = df_ed_raw,
  date_col      = "ENC_ADMIT_DATE",   # or "ENC_ADMIT_DATE"
  community_col = "community", # or "CA_NAME", depending on your data
  label_prefix  = "ED visits"
)

flow_ed

flow_all <- bind_rows(flow_mort, flow_ed)

flow_all

# Write to CSV for documentation / supplement
readr::write_csv(flow_all,
                 "output/tables/flowchart_counts_mortality_ed.csv")



# install.packages("DiagrammeR")  # if needed
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

gr <- grViz("
digraph cohort_flow {

  graph [
    rankdir = LR,
    labelloc = t,
    label = 'Study cohort flow for mortality and ED analyses',
    fontsize = 24,
    fontname = 'Helvetica'
  ]

  node [
    shape = box,
    style = filled,
    color = '#444444',
    fillcolor = '#f7f7f7',
    fontname = 'Helvetica',
    fontsize = 16,
    margin = '0.25,0.20'
  ]

  edge [
    color = '#444444',
    arrowsize = 0.9
  ]

  # ---------- Mortality side ----------
  subgraph cluster_mortality {
    label = 'B. Mortality records';
    fontsize = 20;
    fontname = 'Helvetica';
    color = '#bbbbbb';

    m0 [label = 'Raw death records\\nDeaths_raw\\n n = 653,808'];
    m1 [label = 'Records 2011\u20132022\\n n = 242,626'];
    m2 [label = 'Valid Chicago community area\\n n = 242,626'];
    m3 [label = 'Warm months only\\n(May\u2013September)\\n n = 96,791'];

    m4 [label = 'CVD deaths\\n(any CVD flag = 1)\\n n = 46,843',
        fillcolor = '#e0f3db'];
    m5 [label = 'Non-CVD deaths\\n n = 49,948',
        fillcolor = '#fee0d2'];

    m6 [label = 'CHD deaths\\n n = 15,574',
        fillcolor = '#e0f3db'];
    m7 [label = 'MI deaths\\n n = 5,075',
        fillcolor = '#e0f3db'];
    m8 [label = 'Stroke deaths\\n n = 7,463',
        fillcolor = '#e0f3db'];

    # arrows
    m0 -> m1 -> m2 -> m3;
    m3 -> m4;
    m3 -> m5;
    m4 -> m6;
    m4 -> m7;
    m4 -> m8;
  }

  # ---------- ED side ----------
  subgraph cluster_ed {
    label = 'A. Emergency department visits';
    fontsize = 20;
    fontname = 'Helvetica';
    color = '#bbbbbb';

    e0 [label = 'Raw ED visit records\\n df_ed_raw\\n n = 1,927,640'];
    e1 [label = 'Records 2011\u20132022\\n n = 1,834,968'];
    e2 [label = 'Valid Chicago community area\\n n = 1,834,968'];
    e3 [label = 'Warm months only\\n(May\u2013September)\\n n = 782,416'];

    e4 [label = 'CVD ED visits\\n(any CVD flag = 1)\\n n = 219,044',
        fillcolor = '#e0f3db'];
    e5 [label = 'Non-CVD ED visits\\n n = 563,372',
        fillcolor = '#fee0d2'];

    e6 [label = 'CHD ED visits\\n n = 37,404',
        fillcolor = '#e0f3db'];
    e7 [label = 'MI ED visits\\n n = 4,758',
        fillcolor = '#e0f3db'];
    e8 [label = 'Stroke ED visits\\n n = 19,236',
        fillcolor = '#e0f3db'];

    # arrows
    e0 -> e1 -> e2 -> e3;
    e3 -> e4;
    e3 -> e5;
    e4 -> e6;
    e4 -> e7;
    e4 -> e8;
  }

}
")

# Convert to SVG
svg <- export_svg(gr)

# Save as PNG (high resolution)
rsvg_png(charToRaw(svg),
         file = "flowchart_mortality_ed.png",
         width = 3000,   # adjust size as needed
         height = 2000)

# Save as PDF
rsvg_pdf(charToRaw(svg),
         file = "flowchart_mortality_ed.pdf")



primary_cvd_n <- deaths %>%
  filter(grepl("^I", UNDCAUSE)) %>%
  nrow()

primary_cvd_n





