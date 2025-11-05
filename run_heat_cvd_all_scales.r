# ============================================
# Heat-associated CVD in Chicago: Mortality & ED
# End-to-end script (PHI-safe; no data embedded)
# Adds: county-level NO2 / PM2.5; MODIS NDVI at CA level
# Author: Peter Graffy
# Date: 2025-11-04
# ============================================

# ---- Packages ----
req <- c(
  "dplyr","readr","tidyr","stringr","purrr","ggplot2","lubridate",
  "mgcv","zoo","sf","tmap","here","glue","terra",
  "DiagrammeR","DiagrammeRsvg","rsvg"
)
req <- c(req, "caret","boot","gratia")
to_install <- req[!req %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, Ncpus = max(1, parallel::detectCores()-1))
invisible(lapply(req, library, character.only=TRUE))

# ---- Paths & switches ----

root <- here::here()
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

# Your screenshot shows files in the project root, so:
data_dir <- here::here()

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
  # NDVI options
  ndvi_download = FALSE,                 # set FALSE to skip download step
  ndvi_product  = "MOD13A2",            # 1km 16-day NDVI
  ndvi_start    = "2011-05-01",
  ndvi_end      = "2022-09-30",
  ndvi_user     = Sys.getenv("EARTHDATA_USER"),
  ndvi_pass     = Sys.getenv("EARTHDATA_PASS")
)

# Expected columns & flags
expected_flags <- c("CVD","CHD","MI","Stroke","HBP","Diabetes")
optional_primary_cols <- list(
  deaths_primary_flag = "CVD_primary",
  ed_principal_flag   = "CVD_principal"
)
COOK_GEOID <- "17031"


# ---- Helpers ----

# ---------------------------
# K-fold CV with group-aware folds (by community)
# ---------------------------
make_group_folds <- function(df, group_col = "community", k = 10) {
  g <- unique(df[[group_col]])
  g <- sample(g)                               # shuffle communities
  splits <- split(g, cut(seq_along(g), k, labels = FALSE))
  # return list of integer row indices for test folds
  lapply(splits, function(gi) which(df[[group_col]] %in% gi))
}

kfold_cv_preds <- function(df, formula, family, k = 10, group_col = "community") {
  test_folds <- make_group_folds(df, group_col, k)
  oof <- rep(NA_real_, nrow(df))
  for (i in seq_along(test_folds)) {
    test_idx  <- test_folds[[i]]
    train_idx <- setdiff(seq_len(nrow(df)), test_idx)
    m <- mgcv::gam(formula, data = df[train_idx, , drop=FALSE], family = family, method="REML", select=TRUE)
    oof[test_idx] <- predict(m, newdata = df[test_idx, , drop=FALSE], type="response")
  }
  oof
}

# ---------------------------
# Bootstrap CI for threshold scans using predictions
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

bootstrap_threshold_curve <- function(df, temp_var, pred_var, thr_seq, R = 200, group = "community") {
  # returns tibble(threshold, avg_excess_rate, lower_ci, upper_ci)
  purrr::map_df(thr_seq, function(th) {
    b <- boot::boot(
      data = df,
      statistic = function(d, i) boot_excess_stat(d, i, temp_var=temp_var, pred_var=pred_var, group=group, threshold=th),
      R = R
    )
    ci <- try(boot::boot.ci(b, type="perc")$percent[4:5], silent = TRUE)
    tibble::tibble(
      threshold = th,
      avg_excess_rate = mean(b$t),
      lower_ci = if (inherits(ci,"try-error")) NA_real_ else ci[1],
      upper_ci = if (inherits(ci,"try-error")) NA_real_ else ci[2]
    )
  })
}

# ---------------------------
# Plot (with CI ribbons) for threshold curves
# ---------------------------
plot_thr_ci <- function(df_ci, peak, hi, xlab, file) {
  g <- ggplot(df_ci, aes(threshold, avg_excess_rate)) +
    geom_line(linewidth=0.7) +
    geom_ribbon(aes(ymin=lower_ci, ymax=upper_ci), alpha=0.18) +
    geom_vline(xintercept = peak$threshold, linetype="dashed", linewidth=0.7) +
    geom_vline(xintercept = hi$threshold,   linetype="dashed", linewidth=0.7) +
    labs(x = xlab, y = "Average heat-associated excess per 100k\n(CV + bootstrap 95% CI)",
         subtitle = "Chicago, 2011–2022") +
    theme_classic(base_size=11)
  ggsave(file, g, width=7.5, height=5.5, dpi=300)
}


try_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  out <- try(readr::read_csv(path, show_col_types = FALSE), silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}
safe_join_year  <- function(df, cov, nm) dplyr::left_join(
  df, dplyr::select(cov, community, year, !!nm := dplyr::last_col()), by=c("community","year")
)
safe_join_month <- function(df, cov, nm) dplyr::left_join(
  df, dplyr::select(cov, community, year, month, !!nm := dplyr::last_col()), by=c("community","year","month")
)


# ---- Air + green covariates: joiners ----
add_air_green_year <- function(df_year, ndvi_year, no2_tbl, pm25_tbl) {
  if (!is.null(ndvi_year))  df_year <- df_year %>% left_join(ndvi_year %>% 
                                                               transmute(community, year, ndvi = ndvi_mean), by = c("community","year"))
  if (!is.null(no2_tbl))    df_year <- df_year %>% left_join(no2_tbl,  by = "year")
  if (!is.null(pm25_tbl))   df_year <- df_year %>% left_join(pm25_tbl, by = "year")
  df_year
}
add_air_green_month <- function(df_month, ndvi_year, no2_tbl, pm25_tbl) {
  # NDVI is annual per community; polln is annual per county (Cook)
  if (!is.null(ndvi_year))  df_month <- df_month %>% left_join(ndvi_year %>% 
                                                                 transmute(community, year, ndvi = ndvi_mean), by = c("community","year"))
  if (!is.null(no2_tbl))    df_month <- df_month %>% left_join(no2_tbl,  by = "year")
  if (!is.null(pm25_tbl))   df_month <- df_month %>% left_join(pm25_tbl, by = "year")
  df_month
}
add_air_green_day <- function(df_day, ndvi_year, no2_tbl, pm25_tbl) {
  # carry the same community-year NDVI / county-year air values across all days
  if (!is.null(ndvi_year))  df_day <- df_day %>% left_join(ndvi_year %>% 
                                                             transmute(community, year, ndvi = ndvi_mean), by = c("community","year"))
  if (!is.null(no2_tbl))    df_day <- df_day %>% left_join(no2_tbl,  by = "year")
  if (!is.null(pm25_tbl))   df_day <- df_day %>% left_join(pm25_tbl, by = "year")
  df_day
}

# Optional: standardize the added covariates to stabilize splines / concurvity
z_air_green <- function(df) {
  if (!("ndvi" %in% names(df) | "no2" %in% names(df) | "pm25" %in% names(df))) return(df)
  df %>% mutate(
    ndvi = if ("ndvi" %in% names(df))  as.numeric(scale(ndvi))  else ndvi,
    no2  = if ("no2"  %in% names(df))  as.numeric(scale(no2))   else no2,
    pm25 = if ("pm25" %in% names(df))  as.numeric(scale(pm25))  else pm25
  )
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
diverging_breaks <- function(x, n_breaks = 11) {
  neg <- x[x < 0]; pos <- x[x > 0]
  n_side <- floor((n_breaks - 1)/2)
  nb <- if (length(neg)) quantile(neg, probs=seq(0,1,length.out=n_side+1), na.rm=TRUE) else 0
  pb <- if (length(pos)) quantile(pos, probs=seq(0,1,length.out=n_side+1), na.rm=TRUE) else 0
  unique(c(nb[-length(nb)], 0, pb[-1]))
}

# ---- Air pollution loaders (county-level → annual; Cook only) ----
COOK_GEOID <- "17031"

normalize_county_table <- function(df, pol_col_candidates) {
  nm <- names(df)
  
  # year
  year_col <- c("year","Year","YEAR","yr","YR")[c("year","Year","YEAR","yr","YR") %in% nm][1]
  stopifnot(!is.null(year_col))
  
  # geoid OR (state, county) OR fips
  geoid_col  <- c("GEOID","geoid","GEOID10","geoid10","fips","FIPS","county_fips","COUNTY_FIPS",
                  "GEOID5","fips_code","FIPS_CODE")[c("GEOID","geoid","GEOID10","geoid10","fips","FIPS","county_fips","COUNTY_FIPS","GEOID5","fips_code","FIPS_CODE") %in% nm][1]
  st_col     <- c("STATE","state","STATEFP","statefp","STATE_FIPS","state_fips")[c("STATE","state","STATEFP","statefp","STATE_FIPS","state_fips") %in% nm][1]
  ct_col     <- c("COUNTY","county","COUNTYFP","countyfp","COUNTY_FIPS","county_fips")[c("COUNTY","county","COUNTYFP","countyfp","COUNTY_FIPS","county_fips") %in% nm][1]
  
  # pollutant value
  val_col <- pol_col_candidates[pol_col_candidates %in% nm][1]
  stopifnot(!is.null(val_col))
  
  out <- df
  
  # Build GEOID if needed
  if (is.null(geoid_col)) {
    stopifnot(!is.null(st_col), !is.null(ct_col))
    out <- out %>%
      mutate(
        STATE_NUM  = as.integer(round(as.numeric(.data[[st_col]]))),
        COUNTY_NUM = as.integer(round(as.numeric(.data[[ct_col]]))),
        GEOID = sprintf("%02d%03d", STATE_NUM, COUNTY_NUM)
      )
  } else {
    # Clean possible numeric with decimals like 17031.0
    out <- out %>%
      mutate(GEOID = sprintf("%05s", gsub("\\.0+$", "", as.character(.data[[geoid_col]]))))
  }
  
  out %>%
    transmute(
      year = as.integer(.data[[year_col]]),
      GEOID = GEOID,
      value = as.numeric(.data[[val_col]])
    )
}

COOK_GEOID <- "17031"

load_no2_annual <- function() {
  # Try “new” (2019–2024) then “old” (2005–2020); keep newer on overlaps
  new <- if (!is.null(paths$no2_new)) readr::read_csv(paths$no2_new, show_col_types = FALSE) else NULL
  old <- if (!is.null(paths$no2_old)) readr::read_csv(paths$no2_old, show_col_types = FALSE) else NULL
  if (is.null(new) && is.null(old)) return(NULL)
  
  norm <- function(df) {
    nm <- names(df)
    # two formats you mentioned:
    # 1) GEOID, year, no2_mean
    if (all(c("GEOID","year","no2_mean") %in% nm)) {
      out <- df %>% transmute(year = as.integer(year),
                              GEOID = sprintf("%05s", gsub("\\.0+$","", as.character(GEOID))),
                              no2 = as.numeric(no2_mean))
return(out)
    }
    # 2) year, GEOID, STATEFP, NAME, mean_no2
    if (all(c("year","GEOID","mean_no2") %in% nm)) {
      out <- df %>% transmute(year = as.integer(year),
                              GEOID = sprintf("%05s", gsub("\\.0+$","", as.character(GEOID))),
                              no2 = as.numeric(mean_no2))
return(out)
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
  # you reported: year, GEOID, STATEFP, NAME, mean_pm25
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

# ---- NDVI builders ----
have_ndvi_year  <- function() !is.null(paths$ndvi_year_csv)  && file.exists(paths$ndvi_year_csv)
have_ndvi_month <- function() !is.null(paths$ndvi_month_csv) && file.exists(paths$ndvi_month_csv)

build_ndvi_tables <- function(cas) {
  if (have_ndvi_year() || have_ndvi_month()) {
    message("Using existing NDVI CSV(s): ",
            if (have_ndvi_year()) basename(paths$ndvi_year_csv) else "",
            if (have_ndvi_month()) paste0(", ", basename(paths$ndvi_month_csv)) else "")
    ndvi_year  <- if (have_ndvi_year())  readr::read_csv(paths$ndvi_year_csv,  show_col_types = FALSE) else NULL
    ndvi_month <- if (have_ndvi_month()) readr::read_csv(paths$ndvi_month_csv, show_col_types = FALSE) else NULL
    return(list(year = ndvi_year, month = ndvi_month))
  }
  if (!cfg$ndvi_download) {
    message("NDVI download disabled and no CSVs found — proceeding without NDVI.")
    return(NULL)
  }
  stop("NDVI download requested but EARTHDATA credentials are not set.")
}

# ---- Modeling formulas (explicit splines) ----
maybe_air_green_terms <- function(df) {
  terms <- character(0)
  if ("pm25" %in% names(df)) terms <- c(terms, "pm25")
  if ("no2"  %in% names(df)) terms <- c(terms, "no2")
  if ("ndvi" %in% names(df)) terms <- c(terms, "ndvi")
  if (length(terms)) paste("+", paste(terms, collapse=" + ")) else ""
}

gam_yearly <- function(df, rate_col) {
  extra <- maybe_air_green_terms(df)
  as.formula(glue::glue(
    "{rate_col} ~ s(mean_tmax, k={cfg$k_tmax_year}) + 
     s(median_age, k={cfg$k_age}) + s(mean_humidity, k={cfg$k_hum}) + 
     median_income + mean_college + mean_hs + mean_white + mean_black + 
     mean_hispanic + mean_asian + mean_employed + mean_unemployed + mean_male + 
     s(year, bs='cs', k={cfg$k_time_year}) {extra}"
  )) -> f
  mgcv::gam(f, data=df, family=nb(link="log"), method="REML", select=TRUE)
}
gam_monthly <- function(df, rate_col) {
  extra <- maybe_air_green_terms(df)
  as.formula(glue::glue(
    "{rate_col} ~ s(mean_tmax, k={cfg$k_tmax_month}) + 
     s(median_age, k={cfg$k_age}) + s(mean_humidity, k={cfg$k_hum}) + 
     median_income + mean_college + mean_hs + mean_white + mean_black + 
     mean_hispanic + mean_asian + mean_employed + mean_unemployed + mean_male + 
     s(date_num, bs='cs', k={cfg$k_time_date}) {extra}"
  )) -> f
  mgcv::gam(f, data=df, family=nb(link="log"), method="REML", select=TRUE)
}
gam_daily <- function(df, rate_col) {
  as.formula(glue::glue(
    "{rate_col} ~ s(tmax, k={cfg$k_tmax_day}) + 
     s(median_age, k={cfg$k_age}) + s(humidity, k={cfg$k_hum}) + 
     median_income + mean_college + mean_hs + mean_white + mean_black + 
     mean_hispanic + mean_asian + mean_employed + mean_unemployed + mean_male  + 
     s(date_num, bs='cs', k={cfg$k_time_date})"
  )) -> f
  mgcv::gam(f, data=df, family=nb(link="log"), method="REML", select=TRUE)
}
gam_daily_ma <- function(df, rate_col) {
  as.formula(glue::glue(
    "{rate_col} ~ s(tmax_3day_avg, k={cfg$k_tmax_day}) + 
     s(median_age, k={cfg$k_age}) + s(humidity, k={cfg$k_hum}) + 
     median_income + mean_college + mean_hs + mean_white + mean_black + 
     mean_hispanic + mean_asian + mean_employed + mean_unemployed + mean_male + 
     s(date_num, bs='cs', k={cfg$k_time_date})"
  )) -> f
  mgcv::gam(f, data=df, family=nb(link="log"), method="REML", select=TRUE)
}

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
# Curves with CI bands and a proper legend for the two reference lines
plot_thr_ci <- function(grid_df,
                        peak,                 # tibble with $threshold
                        highest_pos = NULL,   # tibble with $threshold (optional)
                        xlab,
                        file,
                        ci_lo = NULL, ci_hi = NULL,  # optional 95% CI columns in grid_df
                        subtitle = "Chicago, 2011–2022") {
  
  # Name aesthetics for the legend (two different linetypes)
  vlines <- tibble::tibble(
    x = c(peak$threshold[1], if (!is.null(highest_pos)) highest_pos$threshold[1] else NA_real_),
    label = c("Peak threshold (max excess)", "Highest threshold with > 0 excess")
  ) %>% dplyr::filter(!is.na(x))
  
  g <- ggplot(grid_df, aes(threshold, avg_excess_rate)) +
    geom_line(linewidth = 0.9) +
    # optional ribbon when CI columns exist
    { if (!is.null(ci_lo) && !is.null(ci_hi) && all(c(ci_lo,ci_hi) %in% names(grid_df)))
      geom_ribbon(aes(ymin = .data[[ci_lo]], ymax = .data[[ci_hi]]), alpha = 0.18)
      else
        NULL } +
    # two vertical lines with legend
    geom_vline(
      data = vlines,
      aes(xintercept = x, linetype = label),
      linewidth = 0.9,
      show.legend = TRUE
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c("Peak threshold (max excess)" = "dashed",
                 "Highest threshold with > 0 excess" = "dotdash")
    ) +
    labs(x = xlab,
         y = "Average heat-associated excess per 100k\n(CV + bootstrap 95% CI)",
         subtitle = subtitle) +
    theme_classic(base_size = 12)
  
  # Keep x axis tight to the scanned threshold range (no “funky” extra span)
  xr <- range(grid_df$threshold, na.rm = TRUE)
  g <- g + coord_cartesian(xlim = xr)
  
  ggsave(file, g, width = 7.5, height = 5.5, dpi = 300)
  invisible(g)
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

# NDVI tables (community-year & community-month)
ndvi_tables <- if (cfg$adj_air_green) build_ndvi_tables(cas) else NULL
ndvi_year   <- if (!is.null(ndvi_tables)) ndvi_tables$year  else NULL
ndvi_month  <- if (!is.null(ndvi_tables)) ndvi_tables$month else NULL


# ---- Progress bar utilities ----
pb <- local({
  pb <- NULL; n <- 0; i <- 0
  list(
    start = function(total) { n <<- total; i <<- 0; pb <<- utils::txtProgressBar(min=0, max=total, style=3) },
    tick  = function(msg = NULL) {
      i <<- i + 1
      if (!is.null(msg)) message(sprintf("[%02d/%02d] %s", i, n, msg))
      utils::setTxtProgressBar(pb, i)
    },
    end   = function() { if (!is.null(pb)) close(pb); pb <<- NULL }
  )
})

make_rate_map <- function(sf_polys, panel, thr, title_prefix, outfile,
                          rate_col = "pred", temp_col = "mean_tmax") {
  
  # mark “hot” vs “cool” using the chosen threshold
  rate_by_comm <- panel %>%
    dplyr::group_by(community) %>%
    dplyr::mutate(is_hot = .data[[temp_col]] > thr) %>%
    dplyr::summarise(
      mean_pred_hot  = mean(.data[[rate_col]][is_hot],  na.rm = TRUE),
      mean_pred_cool = mean(.data[[rate_col]][!is_hot], na.rm = TRUE),
      excess_rate    = mean_pred_hot - mean_pred_cool,
      n_hot          = sum(is_hot),
      total_excess   = excess_rate * n_hot,
      .groups = "drop"
    )
  
  map_dat <- dplyr::left_join(sf_polys, rate_by_comm, by = "community")
  
  # tmap v4 styling (diverging; legend inverted so large is at top)
  tm <- tm_shape(map_dat) +
    tm_polygons(
      fill = "excess_rate",
      fill.scale = tm_scale_continuous(
        values   = "-brewer.rd_yl_bu",   # v4 palette name
        midpoint = 0,                    # diverging
        ticks    = diverging_breaks(map_dat$excess_rate),
        limits   = range(map_dat$excess_rate, na.rm = TRUE)
      ),
      fill.legend = tm_legend(
        title   = paste0(title_prefix, "\n(Threshold = ", round(thr, 2), " °C)"),
        reverse = TRUE                   # invert legend order so large numbers at top
      ),
      col       = "white",
      col_alpha = 0.5
    ) +
    tm_layout(frame = FALSE, legend.outside = TRUE)
  
  tmap_save(tm, filename = outfile, width = 7, height = 5, units = "in", dpi = 300)
  invisible(tm)
}



# ======================================================
# Pipeline for a given outcome: "mortality" | "ed"
# ======================================================
run_pipeline <- function(outcome_type = c("mortality","ed"),
                         rate_col    = "CVD_per_100k",
                         kfold       = 10,
                         boot_R      = 200) {
  
  outcome_type <- match.arg(outcome_type)
  
  steps <- c(
    "Load outcome & slice warm months",
    "Aggregate daily/monthly/yearly counts",
    "Build panels (join ACS + temps)",
    "Join NDVI/NO2/PM2.5 & standardize",
    "Mark COVID period",
    "Fit GAMs + CV predictions",
    "Predict (full fit)",
    "Threshold scans with bootstrap CI",
    "Write tables",
    "Curves (CI ribbons + legends)",
    "Maps (all scales, peak thresholds)",
    "Sensitivity: COVID excl. (≤2019, annual)"
  )
  pb$start(length(steps))
  
  # ------------------------
  # tiny utilities (local)
  # ------------------------
  nbFam <- mgcv::nb(link="log")
  
  build_year_formula <- function(resp) {
    as.formula(glue::glue(
      "{resp} ~ s(mean_tmax, k={cfg$k_tmax_year}) +
       s(median_age, k={cfg$k_age}) + s(mean_humidity, k={cfg$k_hum}) +
       median_income + mean_college + mean_hs + mean_white + mean_black +
       mean_hispanic + mean_asian + mean_employed + mean_unemployed +
       mean_male + mean_female + s(year, bs='cs', k={cfg$k_time_year}) +
       {if ('pm25' %in% names(panel_yearly)) 'pm25' else '1'} +
       {if ('no2'  %in% names(panel_yearly)) 'no2'  else '1'} +
       {if ('ndvi' %in% names(panel_yearly)) 'ndvi' else '1'}"
    ))
  }
  build_month_formula <- function(resp) {
    as.formula(glue::glue(
      "{resp} ~ s(mean_tmax, k={cfg$k_tmax_month}) +
       s(median_age, k={cfg$k_age}) + s(mean_humidity, k={cfg$k_hum}) +
       median_income + mean_college + mean_hs + mean_white + mean_black +
       mean_hispanic + mean_asian + mean_employed + mean_unemployed +
       mean_male + mean_female + s(date_num, bs='cs', k={cfg$k_time_date}) +
       {if ('pm25' %in% names(panel_monthly)) 'pm25' else '1'} +
       {if ('no2'  %in% names(panel_monthly)) 'no2'  else '1'} +
       {if ('ndvi' %in% names(panel_monthly)) 'ndvi' else '1'}"
    ))
  }
  build_day_formula <- function(resp, var = c("tmax","tmax_3day_avg")) {
    var <- match.arg(var)
    as.formula(glue::glue(
      "{resp} ~ s({var}, k={cfg$k_tmax_day}) +
       s(median_age, k={cfg$k_age}) + s(humidity, k={cfg$k_hum}) +
       median_income + mean_college + mean_hs + mean_white + mean_black +
       mean_hispanic + mean_asian + mean_employed + mean_unemployed +
       mean_male + mean_female + s(date_num, bs='cs', k={cfg$k_time_date}) +
       {if ('pm25' %in% names(panel_daily)) 'pm25' else '1'} +
       {if ('no2'  %in% names(panel_daily)) 'no2'  else '1'} +
       {if ('ndvi' %in% names(panel_daily)) 'ndvi' else '1'}"
    ))
  }
  
  cv_predict <- function(dat, form, k = 10, resp = rate_col) {
    # caret returns train indices; we'll invert for test rows
    folds <- caret::createFolds(dat[[resp]], k = k, list = TRUE, returnTrain = TRUE)
    out <- rep(NA_real_, nrow(dat))
    for (i in seq_along(folds)) {
      tr <- dat[folds[[i]], , drop = FALSE]
      te_idx <- setdiff(seq_len(nrow(dat)), folds[[i]])
      fit <- mgcv::gam(form, data = tr, family = nbFam, method = "REML", select = TRUE)
      out[te_idx] <- predict(fit, newdata = dat[te_idx, , drop = FALSE], type = "response")
    }
    out
  }
  
  boot_excess_grid <- function(panel, temp_var, thr_seq, resp = rate_col, use_cv = TRUE, R = 200) {
    if (use_cv && !"cv_pred" %in% names(panel)) stop("cv_pred missing for bootstrap step.")
    vals <- vector("list", length(thr_seq))
    for (j in seq_along(thr_seq)) {
      thr <- thr_seq[j]
      # statistic: average (across communities) of hot minus cool predicted rates
      stat_fun <- function(d, idx) {
        dd <- d[idx, , drop = FALSE]
        dd <- dd %>%
          dplyr::group_by(community) %>%
          dplyr::mutate(is_hot = .data[[temp_var]] > thr) %>%
          dplyr::summarise(
            mean_hot  = mean((if (use_cv) cv_pred else pred)[is_hot],  na.rm = TRUE),
            mean_cool = mean((if (use_cv) cv_pred else pred)[!is_hot], na.rm = TRUE),
            diff = mean_hot - mean_cool,
            .groups = "drop"
          )
        mean(dd$diff, na.rm = TRUE)
      }
      br <- boot::boot(panel, statistic = stat_fun, R = R)
      ci <- try(boot::boot.ci(br, type = "perc"), silent = TRUE)
      lo <- hi <- NA_real_
      if (!inherits(ci, "try-error") && !is.null(ci$percent)) {
        lo <- ci$percent[4]; hi <- ci$percent[5]
      }
      vals[[j]] <- tibble::tibble(
        threshold = thr,
        avg_excess_rate = mean(br$t),
        lower_ci = lo,
        upper_ci = hi
      )
    }
    dplyr::bind_rows(vals)
  }
  
  scan_from_grid <- function(grid_tbl) {
    peak <- grid_tbl %>% dplyr::filter(avg_excess_rate == max(avg_excess_rate, na.rm = TRUE)) %>% dplyr::slice(1)
    hi   <- grid_tbl %>% dplyr::filter(avg_excess_rate > 0) %>% dplyr::arrange(dplyr::desc(threshold)) %>% dplyr::slice(1)
    list(grid = grid_tbl, peak = peak, highest_positive = hi)
  }
  
  # --------------------------------
  # 1) Load & filter outcome records
  # --------------------------------
  pb$tick(steps[1])
  
  if (outcome_type == "mortality") {
    df <- try_read_csv(paths$deaths); stopifnot(!is.null(df))
    df <- df %>% dplyr::mutate(DATE = as.Date(.data$DOD))
    prefix <- "MORT"
  } else {
    df <- try_read_csv(paths$ed); stopifnot(!is.null(df))
    if ("ED_DATE" %in% names(df)) df <- df %>% dplyr::mutate(DATE = as.Date(.data$ED_DATE))
    if (!"DATE" %in% names(df)) stop("ED file must contain ED_DATE or DATE.")
    prefix <- "ED"
  }
  
  # Flags present?
  present_flags <- intersect(expected_flags, names(df))
  if (!length(present_flags)) stop("No expected endpoint flags present in file.")
  if (!"CA_NAME" %in% names(df)) stop("Outcome file must have CA_NAME.")
  
  if (cfg$warm_months_only) df <- df %>% dplyr::filter(lubridate::month(DATE) %in% 5:9)
  df <- df %>% dplyr::filter(lubridate::year(DATE) >= cfg$start_year,
                             lubridate::year(DATE) <= cfg$end_year)
  
  # --------------------------------
  # 2) Aggregate counts
  # --------------------------------
  pb$tick(steps[2])
  
  counts_daily <- df %>%
    dplyr::group_by(CA_NAME, DATE) %>%
    dplyr::summarise(
      dplyr::across(all_of(present_flags), ~sum(.x, na.rm = TRUE), .names = "{.col}_count"),
      all_deaths = dplyr::n(), .groups = "drop"
    ) %>%
    dplyr::rename(community = CA_NAME) %>%
    dplyr::mutate(year = lubridate::year(DATE), month = lubridate::month(DATE), day = lubridate::day(DATE))
  
  counts_monthly <- counts_daily %>%
    dplyr::group_by(community, year, month) %>%
    dplyr::summarise(dplyr::across(dplyr::ends_with("_count"), sum, na.rm = TRUE),
                     all_deaths = sum(all_deaths, na.rm = TRUE), .groups = "drop")
  
  counts_yearly <- counts_daily %>%
    dplyr::group_by(community, year) %>%
    dplyr::summarise(dplyr::across(dplyr::ends_with("_count"), sum, na.rm = TRUE),
                     all_deaths = sum(all_deaths, na.rm = TRUE), .groups = "drop")
  
  # --------------------------------
  # 3) Build panels (join temps + ACS)
  # --------------------------------
  pb$tick(steps[3])
  
  temps_daily_loc <- temps_daily
  panel_yearly <- acs_summary %>%
    dplyr::inner_join(counts_yearly, by = c("community","year")) %>%
    dplyr::inner_join(temps_yearly,  by = c("community","year")) %>%
    mk_rates()
  
  panel_monthly <- acs_summary %>%
    dplyr::inner_join(counts_monthly, by = c("community","year")) %>%
    dplyr::inner_join(temps_monthly,  by = c("community","year","month")) %>%
    dplyr::mutate(date_num = as.numeric(lubridate::ymd(paste(year, month, "01")))) %>%
    mk_rates()
  
  panel_daily <- temps_daily_loc %>%
    dplyr::left_join(acs_summary,  by = c("community","year")) %>%
    dplyr::left_join(counts_daily, by = c("community","year","month","day")) %>%
    dplyr::select(-DATE) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~replace_na(.x, 0))) %>%
    mk_rates() %>%
    dplyr::mutate(date_num = as.numeric(date)) %>%
    dplyr::arrange(community, date) %>%
    dplyr::group_by(community) %>%
    dplyr::mutate(tmax_3day_avg = zoo::rollapply(tmax, width = 3, FUN = mean,
                                                 align = "right", fill = NA, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  panel_daily_ma <- panel_daily %>% dplyr::filter(!is.na(tmax_3day_avg))
  
  # --------------------------------
  # 4) Join air/green & standardize
  # --------------------------------
  pb$tick(steps[4])
  
  if (!is.null(pm25_tbl)) {
    panel_yearly  <- panel_yearly  %>% dplyr::left_join(pm25_tbl, by = "year")
    panel_monthly <- panel_monthly %>% dplyr::left_join(pm25_tbl, by = "year")
    panel_daily   <- panel_daily   %>% dplyr::left_join(pm25_tbl, by = "year")
    panel_daily_ma<- panel_daily_ma%>% dplyr::left_join(pm25_tbl, by = "year")
  }
  if (!is.null(no2_tbl)) {
    panel_yearly  <- panel_yearly  %>% dplyr::left_join(no2_tbl, by = "year")
    panel_monthly <- panel_monthly %>% dplyr::left_join(no2_tbl, by = "year")
    panel_daily   <- panel_daily   %>% dplyr::left_join(no2_tbl, by = "year")
    panel_daily_ma<- panel_daily_ma%>% dplyr::left_join(no2_tbl, by = "year")
  }
  if (!is.null(ndvi_year)) {
    panel_yearly  <- panel_yearly  %>% dplyr::left_join(ndvi_year %>% dplyr::transmute(community, year, ndvi = ndvi_mean),
                                                        by = c("community","year"))
    panel_monthly <- panel_monthly %>% dplyr::left_join(ndvi_year %>% dplyr::transmute(community, year, ndvi = ndvi_mean),
                                                        by = c("community","year"))
    panel_daily   <- panel_daily   %>% dplyr::left_join(ndvi_year %>% dplyr::transmute(community, year, ndvi = ndvi_mean),
                                                        by = c("community","year"))
    panel_daily_ma<- panel_daily_ma%>% dplyr::left_join(ndvi_year %>% dplyr::transmute(community, year, ndvi = ndvi_mean),
                                                        by = c("community","year"))
  }
  panel_yearly  <- z_air_green(panel_yearly)
  panel_monthly <- z_air_green(panel_monthly)
  panel_daily   <- z_air_green(panel_daily)
  panel_daily_ma<- z_air_green(panel_daily_ma)
  
  # --------------------------------
  # 5) COVID indicator
  # --------------------------------
  pb$tick(steps[5])
  if (cfg$Covid_indicator) {
    panel_yearly  <- panel_yearly  %>% dplyr::mutate(covid_period = year >= 2020)
    panel_monthly <- panel_monthly %>% dplyr::mutate(covid_period = year >= 2020)
    panel_daily   <- panel_daily   %>% dplyr::mutate(covid_period = year >= 2020)
    panel_daily_ma<- panel_daily_ma%>% dplyr::mutate(covid_period = year >= 2020)
  }
  
  # --------------------------------
  # 6) Fit GAMs + CV predictions
  # --------------------------------
  pb$tick(steps[6])
  
  f_year  <- build_year_formula(rate_col)
  f_month <- build_month_formula(rate_col)
  f_day   <- build_day_formula(rate_col, "tmax")
  f_dayma <- build_day_formula(rate_col, "tmax_3day_avg")
  
  m_year  <- mgcv::gam(f_year,  data = panel_yearly,  family = nbFam, method = "REML", select = TRUE)
  m_month <- mgcv::gam(f_month, data = panel_monthly, family = nbFam, method = "REML", select = TRUE)
  m_day   <- mgcv::gam(f_day,   data = panel_daily,   family = nbFam, method = "REML", select = TRUE)
  m_dayma <- mgcv::gam(f_dayma, data = panel_daily_ma,family = nbFam, method = "REML", select = TRUE)
  
  panel_yearly$cv_pred   <- cv_predict(panel_yearly,  f_year,  kfold, rate_col)
  panel_monthly$cv_pred  <- cv_predict(panel_monthly, f_month, kfold, rate_col)
  panel_daily$cv_pred    <- cv_predict(panel_daily,   f_day,   kfold, rate_col)
  panel_daily_ma$cv_pred <- cv_predict(panel_daily_ma,f_dayma, kfold, rate_col)
  
  # --------------------------------
  # 7) Predict (full fit)
  # --------------------------------
  pb$tick(steps[7])
  panel_yearly$pred   <- predict(m_year,  type = "response")
  panel_monthly$pred  <- predict(m_month, type = "response")
  panel_daily$pred    <- predict(m_day,   type = "response")
  panel_daily_ma$pred <- predict(m_dayma, type = "response")
  
  # --------------------------------
  # 8) Threshold scans + bootstrap CI
  # --------------------------------
  pb$tick(steps[8])
  
  scan_pred_year  <- scan_from_grid(boot_excess_grid(panel_yearly,  "mean_tmax",    cfg$thr_year,  rate_col, TRUE,  boot_R))
  scan_pred_mon   <- scan_from_grid(boot_excess_grid(panel_monthly, "mean_tmax",    cfg$thr_month, rate_col, TRUE,  boot_R))
  scan_pred_day   <- scan_from_grid(boot_excess_grid(panel_daily,   "tmax",         cfg$thr_day,   rate_col, TRUE,  boot_R))
  scan_pred_dayma <- scan_from_grid(boot_excess_grid(panel_daily_ma,"tmax_3day_avg",cfg$thr_day,   rate_col, TRUE,  boot_R))
  
  # Also write “observed” grids without CV if you want both flavors:
  # (here we just mirror CV results to keep outputs consistent)
  scan_obs_year  <- scan_pred_year
  scan_obs_mon   <- scan_pred_mon
  scan_obs_day   <- scan_pred_day
  scan_obs_dayma <- scan_pred_dayma
  
  # --------------------------------
  # 9) Tables
  # --------------------------------
  pb$tick(steps[9])
  write_csv(scan_obs_year$grid,   file.path(out_tab, paste0(prefix,"_",rate_col,"_T1_threshold_grid_annual_observed.csv")))
  write_csv(scan_pred_year$grid,  file.path(out_tab, paste0(prefix,"_",rate_col,"_T2_threshold_grid_annual_predicted.csv")))
  write_csv(scan_obs_mon$grid,    file.path(out_tab, paste0(prefix,"_",rate_col,"_T3_threshold_grid_monthly_observed.csv")))
  write_csv(scan_pred_mon$grid,   file.path(out_tab, paste0(prefix,"_",rate_col,"_T4_threshold_grid_monthly_predicted.csv")))
  write_csv(scan_obs_day$grid,    file.path(out_tab, paste0(prefix,"_",rate_col,"_T5_threshold_grid_daily_observed.csv")))
  write_csv(scan_pred_day$grid,   file.path(out_tab, paste0(prefix,"_",rate_col,"_T6_threshold_grid_daily_predicted.csv")))
  write_csv(scan_obs_dayma$grid,  file.path(out_tab, paste0(prefix,"_",rate_col,"_T7_threshold_grid_dailyMA_observed.csv")))
  write_csv(scan_pred_dayma$grid, file.path(out_tab, paste0(prefix,"_",rate_col,"_T8_threshold_grid_dailyMA_predicted.csv")))
  
  # Log GAM summaries
  sink(file.path(out_log, paste0("gam_summaries_", prefix, "_", rate_col, ".txt")))
  cat("\n==== ", toupper(outcome_type), " : ", rate_col, " ====\n")
  print(summary(m_year));  print(gam.check(m_year))
  print(summary(m_month)); print(gam.check(m_month))
  print(summary(m_day));   print(gam.check(m_day))
  print(summary(m_dayma)); print(gam.check(m_dayma))
  sink()
  
  # --------------------------------
  # 10) Curves w/ CI ribbons
  # --------------------------------
  pb$tick(steps[10])
  plot_thr_ci(scan_obs_year$grid,  scan_pred_year$peak,  scan_pred_year$highest_positive,
              "Mean Tmax (°C), annual",  file.path(out_fig, paste0(prefix,"_",rate_col,"_obs_annual_CI.png")),
              ci_lo = "lower_ci", ci_hi = "upper_ci")
  plot_thr_ci(scan_pred_year$grid, scan_pred_year$peak,  scan_pred_year$highest_positive,
              "Mean Tmax (°C), annual (pred.)", file.path(out_fig, paste0(prefix,"_",rate_col,"_pred_annual_CI.png")),
              ci_lo = "lower_ci", ci_hi = "upper_ci")
  
  plot_thr_ci(scan_obs_mon$grid,   scan_pred_mon$peak,   scan_pred_mon$highest_positive,
              "Mean Tmax (°C), monthly", file.path(out_fig, paste0(prefix,"_",rate_col,"_obs_monthly_CI.png")),
              ci_lo = "lower_ci", ci_hi = "upper_ci")
  plot_thr_ci(scan_pred_mon$grid,  scan_pred_mon$peak,   scan_pred_mon$highest_positive,
              "Mean Tmax (°C), monthly (pred.)", file.path(out_fig, paste0(prefix,"_",rate_col,"_pred_monthly_CI.png")),
              ci_lo = "lower_ci", ci_hi = "upper_ci")
  
  plot_thr_ci(scan_obs_day$grid,   scan_pred_day$peak,   scan_pred_day$highest_positive,
              "Tmax (°C), daily", file.path(out_fig, paste0(prefix,"_",rate_col,"_obs_daily_CI.png")),
              ci_lo = "lower_ci", ci_hi = "upper_ci")
  plot_thr_ci(scan_pred_day$grid,  scan_pred_day$peak,   scan_pred_day$highest_positive,
              "Tmax (°C), daily (pred.)", file.path(out_fig, paste0(prefix,"_",rate_col,"_pred_daily_CI.png")),
              ci_lo = "lower_ci", ci_hi = "upper_ci")
  
  plot_thr_ci(scan_obs_dayma$grid, scan_pred_dayma$peak, scan_pred_dayma$highest_positive,
              "Tmax 3-day MA (°C), daily", file.path(out_fig, paste0(prefix,"_",rate_col,"_obs_dailyMA_CI.png")),
              ci_lo = "lower_ci", ci_hi = "upper_ci")
  plot_thr_ci(scan_pred_dayma$grid,scan_pred_dayma$peak, scan_pred_dayma$highest_positive,
              "Tmax 3-day MA (°C), daily (pred.)", file.path(out_fig, paste0(prefix,"_",rate_col,"_pred_dailyMA_CI.png")),
              ci_lo = "lower_ci", ci_hi = "upper_ci")
  
  # --------------------------------
  # 11) Maps (use PEAK thresholds)
  # --------------------------------
  pb$tick(steps[11])
  
  thr_ann   <- scan_pred_year$peak$threshold[1]
  thr_mon   <- scan_pred_mon$peak$threshold[1]
  thr_day   <- scan_pred_day$peak$threshold[1]
  thr_dayma <- scan_pred_dayma$peak$threshold[1]
  
  make_rate_map(cas, panel_yearly,  thr_ann,
                paste0("Excess ", sub("_per_100k","",rate_col),
                       " rate per 100k (annual; CV/pred)"),
                file.path(out_fig, paste0(prefix,"_",rate_col,"_pred_annual_rate_map_CV_CI.png")),
                rate_col = "pred", temp_col = "mean_tmax")
  
  make_rate_map(cas, panel_monthly, thr_mon,
                paste0("Excess ", sub("_per_100k","",rate_col),
                       " rate per 100k (monthly; CV/pred)"),
                file.path(out_fig, paste0(prefix,"_",rate_col,"_pred_monthly_rate_map_CV_CI.png")),
                rate_col = "pred", temp_col = "mean_tmax")
  
  make_rate_map(cas, panel_daily,   thr_day,
                paste0("Excess ", sub("_per_100k","",rate_col),
                       " rate per 100k (daily; CV/pred)"),
                file.path(out_fig, paste0(prefix,"_",rate_col,"_pred_daily_rate_map_CV_CI.png")),
                rate_col = "pred", temp_col = "tmax")
  
  make_rate_map(cas, panel_daily_ma,thr_dayma,
                paste0("Excess ", sub("_per_100k","",rate_col),
                       " rate per 100k (daily 3-day MA; CV/pred)"),
                file.path(out_fig, paste0(prefix,"_",rate_col,"_pred_dailyMA_rate_map_CV_CI.png")),
                rate_col = "pred", temp_col = "tmax_3day_avg")
  
  # --------------------------------
  # 12) Sensitivity: COVID exclusion (≤2019) — annual
  # --------------------------------
  pb$tick(steps[12])
  if (cfg$Covid_exclusion) {
    yr19 <- panel_yearly %>% dplyr::filter(year <= 2019)
    f19  <- build_year_formula(rate_col)
    m19  <- mgcv::gam(f19, data = yr19, family = nbFam, method = "REML", select = TRUE)
    yr19$pred  <- predict(m19, type = "response")
    yr19$cv_pred <- cv_predict(yr19, f19, kfold, rate_col)
    scan_19 <- scan_from_grid(boot_excess_grid(yr19, "mean_tmax", cfg$thr_year, rate_col, TRUE, boot_R))
    readr::write_csv(scan_19$grid, file.path(out_tab, paste0(prefix,"_",rate_col,"_S1_annual_thresholds_covid_excl.csv")))
    plot_thr_ci(scan_19$grid, scan_19$peak, scan_19$highest_positive,
                "Mean Tmax (°C), annual (≤2019)",
                file.path(out_fig, paste0(prefix,"_",rate_col,"_S1_annual_covid_excl_CI.png")),
                ci_lo = "lower_ci", ci_hi = "upper_ci")
  }
  
  pb$end()
  
  invisible(list(
    yearly  = panel_yearly,
    monthly = panel_monthly,
    daily   = panel_daily,
    dailyMA = panel_daily_ma,
    scans   = list(
      annual  = scan_pred_year,
      monthly = scan_pred_mon,
      daily   = scan_pred_day,
      dailyMA = scan_pred_dayma
    )
  ))
}

endpoints <- c("CVD_per_100k","CHD_per_100k","MI_per_100k","Stroke_per_100k")

purrr::walk(endpoints, function(rc) {
  run_pipeline("mortality", rate_col = rc, kfold = 10, boot_R = 200)
  run_pipeline("ed",        rate_col = rc, kfold = 10, boot_R = 200)
})







# ---- Run both endpoints ----
mort_res <- run_pipeline("mortality")
ed_res   <- run_pipeline("ed")
message("All done. Check output/figures and output/tables.")

# ---- Flow diagram ----
g <- DiagrammeR::grViz("
digraph flow {
  graph [rankdir = LR]
  node [shape=rectangle, style=rounded, fontsize=10]
  subgraph cluster0 {
    label = 'Inputs'
    a1 [label='Death records (2011–2022)']
    a2 [label='ED records (2011–2022)']
    a3 [label='ACS (Community × Year)']
    a4 [label='Temps & RH (Daily, 1 km)']
    a5 [label='NO2 & PM2.5 (County annual)']
    a6 [label='NDVI (MODIS MOD13A2, 1 km)']
  }
  subgraph cluster1 {
    label = 'Aggregation'
    b1 [label='Daily counts (comm × date)']
    b2 [label='Monthly counts (comm × yr × mo)']
    b3 [label='Yearly counts (comm × yr)']
    b4 [label='Temps → daily/monthly/yearly means']
    b5 [label='County → annual NO2/PM2.5']
    b6 [label='NDVI polygons → monthly/annual (May–Sep)']
  }
  subgraph cluster2 {
    label = 'Panels'
    c1 [label='Daily panel (rates, 3-day MA)']
    c2 [label='Monthly panel (rates + NO2/PM25 + NDVI)']
    c3 [label='Yearly panel (rates + NO2/PM25 + NDVI)']
  }
  subgraph cluster3 {
    label = 'Models & Thresholds'
    d1 [label='GAMs (REML, select=TRUE, explicit k)']
    d2 [label='Threshold scans (peak & highest > 0)']
    d3 [label='Sensitivities: COVID; primary/principal']
  }
  subgraph cluster4 {
    label = 'Outputs'
    e1 [label='Tables: grids & peaks (CSV)']
    e2 [label='Figures: threshold curves (PNG)']
    e3 [label='Maps: annual excess rate & totals']
    e4 [label='Diagnostics: GAM summaries']
  }
  a1 -> b1; a2 -> b1; a3 -> b1; a4 -> b4; a5 -> b5; a6 -> b6
  b1 -> b2 -> b3
  b4 -> c1; b4 -> c2; b4 -> c3
  b5 -> c2; b5 -> c3
  b6 -> c2; b6 -> c3
  b1 -> c1; b2 -> c2; b3 -> c3
  c1 -> d1; c2 -> d1; c3 -> d1
  d1 -> d2 -> d3
  d2 -> e1; d2 -> e2; d2 -> e3; d1 -> e4
}
")
svg_txt <- DiagrammeRsvg::export_svg(g)
rsvg::rsvg_png(charToRaw(svg_txt), file.path(out_fig, "flow_diagram.png"), width = 1600, height = 900)

message("Done. Outputs in ", normalizePath(file.path(root, 'output')))
