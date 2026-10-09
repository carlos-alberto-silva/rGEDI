check_url <- function(url) {
  # Try HEAD without following redirects
  req <- httr2::request(url)
  req <- httr2::req_method(req, "HEAD")
  req <- httr2::req_timeout(req, 10)
  req <- httr2::req_options(req, followlocation = 0)
  req <- httr2::req_user_agent(req, "R (testthat url-check)")

  resp <- try(httr2::req_perform(req), silent = TRUE)

  # If HEAD completely failed, try GET
  if (inherits(resp, "try-error") ||
        httr2::resp_status(resp) %in% c(403, 405)) {
    req <- httr2::request(url)
    req <- httr2::req_method(req, "GET")
    req <- httr2::req_timeout(req, 10)
    req <- httr2::req_options(req, followlocation = 0)
    req <- httr2::req_user_agent(req, "R (testthat url-check)")

    resp <- try(httr2::req_perform(req), silent = TRUE)
  }

  resp
}

test_that("reachable URLs in package sources do not redirect", {
  skip_on_cran()       # CRAN machines don’t like network traffic
  skip_if_offline()    # testthat helper, skips if no internet

  # Use package root whether installed or source
  pkg_root <- testthat::test_path("..", "..")

  # Collect all R/ source files
  r_files <- list.files(file.path(pkg_root),
                        pattern = "\\.R|README\\.md|DESCRIPTION$",
                        full.names = TRUE,
                        recursive = TRUE)

  if (length(r_files) == 0) {
    skip("No .R files found to scan for URLs.")
  }

  # Regex for URLs (basic but avoids trailing punctuation)
  url_pattern <- "(?<!'|\")https?://[^[:space:]'\"<>)}]+"

  urls_found <- unique(unlist(lapply(r_files, function(file) {
    content <- readLines(file, warn = FALSE)
    stringr::str_extract_all(content, url_pattern) |> unlist()
  })))

  # Clean up trailing punctuation if any slipped through
  urls_found <- gsub("[)}.,]+$", "", urls_found)


  if (length(urls_found) == 0) {
    skip("No URLs found in .R files.")
  }

  skip_domains <- c("stackexchange.com")
  for (url in urls_found) {
    if (any(grepl(paste(skip_domains, collapse="|"), url, ignore.case=TRUE))) {
      succeed(sprintf("Skipped URL check for %s (domain blocks bots).", url))
      next
    }
    # Each URL gets its own expectation, not its own test_that()
    resp <- check_url(url)

    if (inherits(resp, "try-error")) {
      succeed(sprintf("URL '%s' could not be reached in this environment.", url))
      next
    }

    status_code <- httr2::resp_status(resp)

    if (status_code >= 300 && status_code < 400) {
      if (startsWith(url,'https://doi.org')) {
        succeed(sprintf("URL '%s' is a DOI URL (status %d).", url, status_code))
        next
      }
      if (startsWith(tolower(url),'https://cran.r-project.org/package=')) {
        succeed(sprintf("URL '%s' is a CRAN package URL (status %d).", url, status_code))
        next
      }
      redirect_url <- httr2::resp_header(resp, "location")
      if (is.null(redirect_url)) redirect_url <- "unknown"
      fail(sprintf("URL '%s' redirects to '%s' (status %d).",
                   url, redirect_url, status_code))
    } else if (status_code %in% c(403, 429, 500, 502, 503, 504)) {
      succeed(sprintf("URL '%s' returned transient status %d.", url, status_code))
    } else if (status_code >= 400) {
      fail(sprintf("URL '%s' returned error status %d.", url, status_code))
    } else {
      succeed(sprintf("URL '%s' is valid (status %d).", url, status_code))
    }
  }
})
