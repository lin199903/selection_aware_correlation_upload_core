options(stringsAsFactors = FALSE)

project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    script <- sub("^--file=", "", file_arg)
    return(dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE))))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT <- project_root()

p <- function(...) file.path(ROOT, ...)

required_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         ". Install them outside the analysis scripts.")
  }
  invisible(TRUE)
}

read_spec <- function() {
  required_packages("jsonlite")
  jsonlite::read_json(p("config", "analysis_spec.json"), simplifyVector = TRUE)
}

ensure_dirs <- function(...) {
  dirs <- unlist(list(...), use.names = FALSE)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

strip_field <- function(x, field) {
  trimws(sub(paste0("^", field, "\\s*:\\s*"), "", as.character(x), ignore.case = TRUE))
}

write_csv_atomic <- function(x, path, row.names = FALSE) {
  tmp <- paste0(path, ".tmp")
  utils::write.csv(x, tmp, row.names = row.names, na = "")
  if (!file.rename(tmp, path)) stop("Could not atomically write ", path)
  invisible(path)
}

write_json_atomic <- function(x, path) {
  required_packages("jsonlite")
  tmp <- paste0(path, ".tmp")
  jsonlite::write_json(x, tmp, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
  if (!file.rename(tmp, path)) stop("Could not atomically write ", path)
  invisible(path)
}

sha256_file <- function(path) {
  required_packages("digest")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

save_session_info <- function(name) {
  ensure_dirs(p("logs"))
  writeLines(capture.output(sessionInfo()), p("logs", paste0(name, "_session_info.txt")))
}

collapse_counts_by_symbol <- function(counts, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols)
  counts <- counts[keep, , drop = FALSE]
  symbols <- unname(symbols[keep])
  out <- rowsum(counts, group = symbols, reorder = FALSE, na.rm = TRUE)
  storage.mode(out) <- "numeric"
  out
}

finite_or_na <- function(x) {
  x[!is.finite(x)] <- NA_real_
  x
}

