

# =========================================================
# NOAA NCEI NDVI (AVHRR-Land) bulk downloader
# Years: 2011–2024
# Source index: https://www.ncei.noaa.gov/data/land-normalized-difference-vegetation-index/access/
# Output: data/ndvi/NCEI/<year>/*.nc + manifest CSV
# PHI-safe
# =========================================================

# ---- Packages ----
req <- c("httr", "rvest", "xml2", "readr", "dplyr", "stringr", "fs", "tools")
to_install <- req[!req %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, Ncpus = max(1, parallel::detectCores()-1))
invisible(lapply(req, library, character.only = TRUE))

# ---- Config ----
base_url <- "https://www.ncei.noaa.gov/data/land-normalized-difference-vegetation-index/access"
years    <- 2011:2024
out_dir  <- fs::path("data", "ndvi", "NCEI")
fs::dir_create(out_dir, recurse = TRUE)

# Download settings
max_retries <- 5
sleep_between <- 0.25      # seconds between requests to be polite
timeout_sec <- 120         # per download attempt

# ---- Helpers ----

# Get directory listing for a given year and return absolute URLs to .nc files
list_year_files <- function(year) {
  idx_url <- sprintf("%s/%d/", base_url, year)
  pg <- try(httr::GET(idx_url, httr::timeout(timeout_sec)), silent = TRUE)
  if (inherits(pg, "try-error") || httr::http_error(pg)) {
    warning("Failed to read index for year: ", year)
    return(character(0))
  }
  # Parse links
  html <- xml2::read_html(httr::content(pg, "raw"))
  links <- rvest::html_attr(rvest::html_nodes(html, "a"), "href")
  # Keep .nc files
  files <- links[grepl("\\.nc$", links, ignore.case = TRUE)]
  if (!length(files)) return(character(0))
  file_urls <- paste0(idx_url, files)
  return(file_urls)
}

# HEAD request to obtain content-length (returns NA if not available)
get_remote_size <- function(url) {
  h <- try(httr::HEAD(url, httr::timeout(timeout_sec)), silent = TRUE)
  if (inherits(h, "try-error") || httr::http_error(h)) return(NA_real_)
  cl <- httr::headers(h)[["content-length"]]
  if (is.null(cl)) return(NA_real_)
  as.numeric(cl)
}

# Download with retries; returns TRUE on success
download_with_retries <- function(url, dest_file) {
  for (i in seq_len(max_retries)) {
    try({
      r <- httr::GET(url,
                     httr::write_disk(dest_file, overwrite = TRUE),
                     httr::progress(),
                     httr::timeout(timeout_sec))
      if (!httr::http_error(r)) {
        return(TRUE)
      }
    }, silent = TRUE)
    # backoff
    Sys.sleep(sleep_between * i)
  }
  FALSE
}

# Ensure directory exists for a given year
ensure_year_dir <- function(year) {
  ydir <- fs::path(out_dir, as.character(year))
  fs::dir_create(ydir)
  ydir
}

# ---- Main ----
results <- list()

for (yr in years) {
  message("Year ", yr, " — listing…")
  urls <- list_year_files(yr)
  if (!length(urls)) {
    warning("No files found for year ", yr, " (skipping).")
    next
  }
  
  ydir <- ensure_year_dir(yr)
  
  for (u in urls) {
    # local file name is basename of URL
    fname <- basename(u)
    dest  <- fs::path(ydir, fname)
    
    # Remote size (may be NA if server doesn't provide)
    remote_size <- get_remote_size(u)
    
    # Skip if already downloaded with same size
    if (fs::file_exists(dest) && !is.na(remote_size)) {
      local_size <- fs::file_size(dest)
      if (as.numeric(local_size) == remote_size) {
        message("✓ Skipping (already present): ", fname)
        results[[length(results) + 1]] <- data.frame(
          year = yr, file = as.character(dest), url = u,
          status = "skipped_exists_size_match",
          size_bytes = as.numeric(local_size),
          md5 = tools::md5sum(dest),
          stringsAsFactors = FALSE
        )
        next
      }
    }
    
    message("↓ Downloading: ", fname)
    ok <- download_with_retries(u, dest)
    Sys.sleep(sleep_between)
    
    if (!ok) {
      warning("Failed to download: ", u)
      results[[length(results) + 1]] <- data.frame(
        year = yr, file = as.character(dest), url = u,
        status = "failed",
        size_bytes = if (fs::file_exists(dest)) as.numeric(fs::file_size(dest)) else NA_real_,
        md5 = if (fs::file_exists(dest)) tools::md5sum(dest) else NA_character_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Validate size if available
    sz <- as.numeric(fs::file_size(dest))
    if (!is.na(remote_size) && !is.na(sz) && remote_size != sz) {
      warning("Size mismatch for ", fname, " (remote=", remote_size, ", local=", sz, "). Keeping file but marking mismatch.")
      status <- "downloaded_size_mismatch"
    } else {
      status <- "downloaded_ok"
    }
    
    results[[length(results) + 1]] <- data.frame(
      year = yr, file = as.character(dest), url = u,
      status = status,
      size_bytes = sz,
      md5 = tools::md5sum(dest),
      stringsAsFactors = FALSE
    )
  }
}

manifest <- if (length(results)) dplyr::bind_rows(results) else {
  tibble::tibble(year=integer(), file=character(), url=character(), status=character(),
                 size_bytes=double(), md5=character())
}

fs::dir_create(fs::path(out_dir, "_manifests"))
readr::write_csv(manifest, fs::path(out_dir, "_manifests", "ndvi_ncei_manifest_2011_2024.csv"))

message("Done. Files are in ", fs::path_abs(out_dir))



