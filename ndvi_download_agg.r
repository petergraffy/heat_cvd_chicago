# =========================================================
# NOAA NCEI NDVI (AVHRR-Land) downloader — July 1 snapshots
# Years: 2011–2024 (configurable)
# Source index: https://www.ncei.noaa.gov/data/land-normalized-difference-vegetation-index/access/
# Output: data/ndvi/NCEI/<year>/<file>.nc + manifest CSV
# =========================================================

# ---- Packages ----
req <- c("httr", "rvest", "xml2", "readr", "dplyr", "stringr", "fs", "tools", "lubridate")
to_install <- req[!req %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, Ncpus = max(1, parallel::detectCores()-1))
invisible(lapply(req, library, character.only = TRUE))

# ---- Config ----
base_url <- "https://www.ncei.noaa.gov/data/land-normalized-difference-vegetation-index/access"
years    <- 2011:2024                    # change if needed
target_md <- "07-01"                     # Month-day string for July 1
out_dir  <- fs::path("data", "ndvi", "NCEI")
fs::dir_create(out_dir, recurse = TRUE)

# polite + robust networking
max_retries   <- 5
sleep_between <- 0.25                    # seconds
timeout_sec   <- 120

# ---- Helpers ----

# Build regex candidates for a given year and date (covers common filename date encodings)
build_date_regex <- function(yr, date) {
  y <- as.integer(yr)
  # YYYY-MM-DD and YYYYMMDD
  iso <- format(date, "%Y-%m-%d")
  compact <- format(date, "%Y%m%d")
  # YYYY-DOY and YYYYDOY (day of year)
  doy <- lubridate::yday(date)
  doy3 <- sprintf("%03d", doy)
  list(
    # exact common patterns
    paste0(y, "-", doy3),            # YYYY-DOY
    paste0(y, doy3),                 # YYYYDOY
    iso,                             # YYYY-MM-DD
    compact                          # YYYYMMDD
  )
}

# Get directory listing for a given year
list_year_links <- function(year) {
  idx_url <- sprintf("%s/%d/", base_url, year)
  pg <- try(httr::GET(idx_url, httr::timeout(timeout_sec)), silent = TRUE)
  if (inherits(pg, "try-error") || httr::http_error(pg)) {
    warning("Failed to read index for year: ", year)
    return(tibble::tibble(href=character(), url=character()))
  }
  html <- xml2::read_html(httr::content(pg, "raw"))
  hrefs <- rvest::html_attr(rvest::html_nodes(html, "a"), "href")
  hrefs <- hrefs[grepl("\\.nc$", hrefs, ignore.case = TRUE)]
  tibble::tibble(
    href = hrefs,
    url  = paste0(idx_url, hrefs)
  )
}

# HEAD → content-length
get_remote_size <- function(url) {
  h <- try(httr::HEAD(url, httr::timeout(timeout_sec)), silent = TRUE)
  if (inherits(h, "try-error") || httr::http_error(h)) return(NA_real_)
  cl <- httr::headers(h)[["content-length"]]
  if (is.null(cl)) return(NA_real_)
  as.numeric(cl)
}

# GET with retries
download_with_retries <- function(url, dest_file) {
  for (i in seq_len(max_retries)) {
    ok <- try({
      r <- httr::GET(url,
                     httr::write_disk(dest_file, overwrite = TRUE),
                     httr::progress(),
                     httr::timeout(timeout_sec))
      !httr::http_error(r)
    }, silent = TRUE)
    if (isTRUE(ok)) return(TRUE)
    Sys.sleep(sleep_between * i)
  }
  FALSE
}

# ---- Main ----
results <- list()

for (yr in years) {
  message("Year ", yr, " — locating July 1 file…")
  links <- list_year_links(yr)
  if (!nrow(links)) {
    warning("No links for year ", yr, "; skipping.")
    next
  }
  
  # Build date candidates for matching
  date_target <- as.Date(paste0(yr, "-", target_md))
  patterns <- build_date_regex(yr, date_target)
  
  # Choose file whose name contains one of the patterns (first match wins)
  idx <- which(vapply(patterns, function(p) any(grepl(p, links$href)), logical(1)))
  chosen <- NULL
  if (length(idx)) {
    # Use the earliest pattern that matches; pick the matching row
    p <- patterns[min(idx)]
    chosen <- links %>% dplyr::filter(grepl(p, href)) %>% slice_head(n = 1)
  } else {
    warning("Could not find July 1 by filename pattern for ", yr, ". Trying fallback: first file in July.")
    # Fallback heuristic: find anything with '-07-' or '07' month in any common encodings
    # (YYYY-07-*, YYYY07*, or DOY around 182/183)
    doy <- lubridate::yday(date_target)
    near_doys <- sprintf("%03d", unique(c(doy-1, doy, doy+1)))  # handle leap years/phase
    fallback_regex <- paste0(
      "(",
      paste(
        c(paste0(yr, "-07-"), paste0(yr, "07"), paste0(yr, "-", near_doys), paste0(yr, near_doys)),
        collapse="|"
      ),
      ")"
    )
    q <- links %>% dplyr::filter(grepl(fallback_regex, href)) %>% arrange(href) %>% slice_head(n=1)
    if (nrow(q)) chosen <- q else {
      warning("No July candidate found for ", yr, ". Skipping year.")
      next
    }
  }
  
  ydir <- fs::path(out_dir, as.character(yr))
  fs::dir_create(ydir)
  fname <- basename(chosen$url)
  dest  <- fs::path(ydir, fname)
  
  # Skip identical-size files
  remote_size <- get_remote_size(chosen$url)
  if (fs::file_exists(dest) && !is.na(remote_size)) {
    local_size <- as.numeric(fs::file_size(dest))
    if (local_size == remote_size) {
      message("✓ Skipping (present, size match): ", fname)
      results[[length(results)+1]] <- data.frame(
        year = yr, file = as.character(dest), url = chosen$url,
        status = "skipped_exists_size_match",
        size_bytes = local_size, md5 = tools::md5sum(dest),
        stringsAsFactors = FALSE
      )
      next
    }
  }
  
  message("↓ Downloading: ", fname)
  ok <- download_with_retries(chosen$url, dest)
  Sys.sleep(sleep_between)
  
  status <- if (ok) "downloaded_ok" else "failed"
  sz <- if (fs::file_exists(dest)) as.numeric(fs::file_size(dest)) else NA_real_
  
  # If we have a remote size, flag mismatches
  if (ok && !is.na(remote_size) && !is.na(sz) && remote_size != sz) {
    warning("Size mismatch for ", fname, " (remote=", remote_size, ", local=", sz, ").")
    status <- "downloaded_size_mismatch"
  }
  
  results[[length(results)+1]] <- data.frame(
    year = yr, file = as.character(dest), url = chosen$url,
    status = status, size_bytes = sz,
    md5 = if (fs::file_exists(dest)) tools::md5sum(dest) else NA_character_,
    stringsAsFactors = FALSE
  )
}

manifest <- if (length(results)) dplyr::bind_rows(results) else {
  tibble::tibble(year=integer(), file=character(), url=character(), status=character(),
                 size_bytes=double(), md5=character())
}

fs::dir_create(fs::path(out_dir, "_manifests"))
readr::write_csv(manifest, fs::path(out_dir, "_manifests", paste0("ndvi_ncei_manifest_", min(years), "_", max(years), "_july1.csv")))

message("Done. Files in: ", fs::path_abs(out_dir))


# =========================================================
# Aggregate NOAA NCEI NDVI (July 1 snapshots) to Chicago
# community areas (CA) and export a CSV for modeling.
#
# Expects files like:
#   data/ndvi/NCEI/2011/AVHRR-Land_v005_AVH13C1_NOAA-19_20110701_*.nc
#
# Output:
#   data/ndvi/derived/ndvi_year_ca.csv  (community, year, ndvi_mean, ndvi_sd, n_cells)
# =========================================================

# ---- Packages ----
req <- c("terra","sf","dplyr","readr","stringr","fs")
to_install <- req[!req %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, Ncpus = max(1, parallel::detectCores()-1))
invisible(lapply(req, library, character.only=TRUE))

# ---- Paths ----
ndvi_root <- fs::path("data","ndvi","NCEI")
cas_path  <- fs::path("data","comm_areas.geojson")   # change if your file lives elsewhere
out_dir   <- fs::path("data","ndvi","derived"); fs::dir_create(out_dir, recurse=TRUE)
out_csv   <- fs::path(out_dir, "ndvi_year_ca.csv")

# ---- Load Chicago community areas ----
stopifnot(fs::file_exists(cas_path))
cas <- sf::st_read(cas_path, quiet = TRUE)
# Expect a "community" column; if not, try common alternatives
if (!"community" %in% names(cas)) {
  alt <- intersect(names(cas), c("community_area","COMMUNITY","ca_name","CA_NAME","COMMUNITY_A"))
  if (length(alt)) cas <- dplyr::rename(cas, community = !!alt[1]) else
    stop("Could not find a 'community' column in the community areas file.")
}
cas <- sf::st_make_valid(cas)

# ---- Find all July-1 NDVI files by year ----
if (!fs::dir_exists(ndvi_root)) stop("NDVI root directory not found: ", ndvi_root)
year_dirs <- fs::dir_ls(ndvi_root, type = "directory", recurse = FALSE)
if (!length(year_dirs)) stop("No yearly subfolders under: ", ndvi_root)

# Helper: pick the NDVI layer from a NetCDF (some files may have multiple layers)
pick_ndvi_layer <- function(r) {
  nms <- names(r)
  # Prefer names that contain 'ndvi' (case-insensitive)
  idx <- which(grepl("ndvi", nms, ignore.case = TRUE))
  if (length(idx)) return(r[[idx[1]]])
  # Otherwise, if there is only one layer, use it; else just take the first
  if (terra::nlyr(r) == 1) r[[1]] else r[[1]]
}

# Process each year → one NDVI mean per CA
all_rows <- list()

for (yd in sort(year_dirs)) {
  yr <- basename(yd)
  if (!grepl("^\\d{4}$", yr)) next
  files <- fs::dir_ls(yd, glob = "*.nc", recurse = FALSE)
  if (!length(files)) {
    message("No .nc files in ", yd, " — skipping.")
    next
  }
  
  # Prefer a file that clearly encodes July 1 in filename; fall back to first .nc
  # Matches: YYYYMMDD, YYYY-MM-DD, YYYYDOY, YYYY-DOY around July 1
  date_target <- as.Date(paste0(yr,"-07-01"))
  cand_pat <- c(
    format(date_target, "%Y%m%d"),
    format(date_target, "%Y-%m-%d"),
    paste0(yr, sprintf("%03d", as.integer(strftime(date_target, "%j")))),
    paste0(yr, "-", sprintf("%03d", as.integer(strftime(date_target, "%j"))))
  )
  pick <- NULL
  for (p in cand_pat) {
    hit <- files[grepl(p, basename(files))]
    if (length(hit)) { pick <- hit[1]; break }
  }
  if (is.null(pick)) {
    # try broader July heuristic
    doy <- as.integer(strftime(date_target, "%j"))
    near <- sprintf("%03d", unique(c(doy-1,doy,doy+1)))
    july_pat <- paste0("(", paste(c(
      paste0(yr, "-07-"), paste0(yr, "07"),
      paste0(yr, "-", near), paste0(yr, near)
    ), collapse="|"), ")")
    hit <- files[grepl(july_pat, basename(files))]
    pick <- if (length(hit)) hit[1] else files[1]
  }
  
  message("Year ", yr, " → using file: ", basename(pick))
  
  # Read raster; auto-scale typically applied if metadata has scale/offset
  r_all <- try(terra::rast(pick), silent = TRUE)
  if (inherits(r_all, "try-error")) {
    warning("Failed to read raster for ", pick, "; skipping year ", yr)
    next
  }
  r <- pick_ndvi_layer(r_all)
  
  # Align CRS & extent
  cas_proj <- sf::st_transform(cas, terra::crs(r))
  cas_v    <- terra::vect(cas_proj)
  
  # Extract per-CA NDVI (mean & sd), removing out-of-range values
  # Clamp to plausible NDVI range [-1, 1] and set others to NA
  r[r < -1 | r > 1] <- NA
  
  ext_mean <- terra::extract(r, cas_v, fun = mean, na.rm = TRUE)
  ext_sd   <- terra::extract(r, cas_v, fun = sd,   na.rm = TRUE)
  ext_n    <- terra::extract(!is.na(r), cas_v, fun = sum, na.rm = TRUE)  # number of contributing cells
  
  # Assemble tidy table
  df <- dplyr::tibble(
    community = cas$community,
    year      = as.integer(yr),
    ndvi_mean = as.numeric(ext_mean[[2]]),
    ndvi_sd   = as.numeric(ext_sd[[2]]),
    n_cells   = as.integer(ext_n[[2]])
  )
  
  all_rows[[length(all_rows) + 1]] <- df
}

ndvi_year_ca <- if (length(all_rows)) dplyr::bind_rows(all_rows) else {
  stop("No NDVI results were produced. Check file locations and formats.")
}

# Optional: basic QC filters / warnings
bad <- ndvi_year_ca %>% dplyr::filter(is.na(ndvi_mean) | n_cells == 0)
if (nrow(bad)) {
  warning("Some community-year rows had no valid NDVI pixels (n=", nrow(bad), "). They will remain as NA.")
}

# Write CSV
readr::write_csv(ndvi_year_ca, out_csv)
message("Wrote: ", fs::path_abs(out_csv))

# ---- Handle missing NDVI (temporal interpolation) ----
message("Interpolating missing NDVI values by community...")

ndvi_interp <- ndvi_year_ca %>%
  dplyr::group_by(community) %>%
  dplyr::arrange(year, .by_group = TRUE) %>%
  dplyr::mutate(
    ndvi_mean_interp = {
      x <- year
      y <- ndvi_mean
      # Identify non-missing years
      ok <- !is.na(y)
      if (sum(ok) >= 2) {
        # Linear interpolation between known values
        approx(x[ok], y[ok], xout = x, rule = 2)$y  # rule=2 allows extrapolation at ends
      } else if (sum(ok) == 1) {
        rep(y[ok], length(y))  # single known value → repeat across all years
      } else {
        rep(NA_real_, length(y))  # no data → all NA
      }
    }
  ) %>%
  dplyr::ungroup()

# Replace ndvi_mean with interpolated version
ndvi_final <- ndvi_interp %>%
  dplyr::mutate(ndvi_mean = ndvi_mean_interp) %>%
  dplyr::select(-ndvi_mean_interp)

# ---- Save final CSV ----
out_csv_final <- fs::path(out_dir, "ndvi_year_ca_interpolated.csv")
readr::write_csv(ndvi_final, out_csv_final)

message("✅ Interpolated NDVI file saved:")
message("   ", fs::path_abs(out_csv_final))

# Quick summary of fill
summary_tbl <- ndvi_final %>%
  summarise(
    n_total = n(),
    n_filled = sum(is.na(ndvi_mean) & !is.na(ndvi_mean)), # sanity check (should be 0 now)
    n_original_na = sum(is.na(ndvi_mean)), 
    ndvi_min = min(ndvi_mean, na.rm=TRUE),
    ndvi_max = max(ndvi_mean, na.rm=TRUE)
  )
print(summary_tbl)












