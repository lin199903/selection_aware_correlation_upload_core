script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("GEOquery", "Biobase", "data.table", "jsonlite", "digest"))

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(data.table)
})

ensure_dirs(p("data", "derived"), p("logs"))

raw_files <- c(
  "GSE126848_Gene_counts_raw.txt.gz",
  "GSE126848_series_matrix.txt.gz",
  "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz",
  "GSE130970_series_matrix.txt.gz",
  "GSE135917_series_matrix.txt.gz",
  "GPL6244.soft.gz"
)
raw_paths <- p("data", "raw", raw_files)
if (!all(file.exists(raw_paths))) {
  stop("Missing raw inputs: ", paste(raw_files[!file.exists(raw_paths)], collapse = ", "))
}
if (any(grepl("_TPM_", basename(raw_paths), ignore.case = TRUE))) {
  stop("TPM input is prohibited in the authoritative raw manifest")
}

manifest <- data.frame(
  file = raw_files,
  bytes = unname(file.info(raw_paths)$size),
  sha256 = vapply(raw_paths, sha256_file, character(1)),
  source = c(
    "NCBI GEO GSE126848 supplementary file",
    "NCBI GEO GSE126848 series matrix",
    "NCBI GEO GSE130970 supplementary tximport counts file",
    "NCBI GEO GSE130970 series matrix",
    "NCBI GEO GSE135917 series matrix",
    "NCBI GEO GPL6244 platform SOFT"
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(manifest, p("data", "derived", "raw_input_manifest.csv"))

read_geo_metadata <- function(file) {
  eset <- getGEO(filename = file, getGPL = FALSE)
  if (is.list(eset)) eset <- eset[[1L]]
  pData(eset)
}

p126 <- read_geo_metadata(p("data", "raw", "GSE126848_series_matrix.txt.gz"))
disease126 <- strip_field(p126[["disease:ch1"]], "disease")
sex126 <- strip_field(p126[["gender:ch1"]], "gender")
group126 <- ifelse(tolower(disease126) %in% c("healthy", "normal", "obese"), "Control",
                   ifelse(tolower(disease126) %in% c("nafld", "nash"), "MASLD", NA_character_))
meta126 <- data.frame(
  sample_id = rownames(p126),
  title = as.character(p126$title),
  raw_column_id = trimws(as.character(p126$description)),
  disease_label = disease126,
  group = group126,
  sex = sex126,
  stringsAsFactors = FALSE
)

raw126_head <- fread(p("data", "raw", "GSE126848_Gene_counts_raw.txt.gz"),
                     header = FALSE, nrows = 1L, showProgress = FALSE)
raw126_dims <- dim(fread(p("data", "raw", "GSE126848_Gene_counts_raw.txt.gz"),
                         header = FALSE, select = 1L, showProgress = FALSE))
stopifnot(ncol(raw126_head) - 1L == nrow(meta126))
raw126_ids <- as.character(unlist(raw126_head[1, -1, with = FALSE], use.names = FALSE))
stopifnot(setequal(raw126_ids, meta126$raw_column_id))
meta126 <- meta126[match(raw126_ids, meta126$raw_column_id), , drop = FALSE]
stopifnot(!anyNA(meta126$sample_id), identical(meta126$raw_column_id, raw126_ids))

p130 <- read_geo_metadata(p("data", "raw", "GSE130970_series_matrix.txt.gz"))
num_field <- function(name, field) as.numeric(strip_field(p130[[name]], field))
steatosis130 <- num_field("steatosis grade:ch1", "steatosis grade")
fibrosis130 <- num_field("fibrosis stage:ch1", "fibrosis stage")
ballooning130 <- num_field("cytological ballooning grade:ch1", "cytological ballooning grade")
inflammation130 <- num_field("lobular inflammation grade:ch1", "lobular inflammation grade")
nas130 <- num_field("nafld activity score:ch1", "nafld activity score")
age130 <- num_field("age at biopsy:ch1", "age at biopsy")
sex130 <- strip_field(p130[["Sex:ch1"]], "Sex")
control130 <- steatosis130 == 0 & fibrosis130 == 0 & ballooning130 == 0
group130 <- ifelse(control130, "Control", "MASLD")
meta130 <- data.frame(
  sample_id = rownames(p130),
  title = as.character(p130$title),
  group = group130,
  sex = sex130,
  age = age130,
  steatosis_grade = steatosis130,
  lobular_inflammation_grade = inflammation130,
  ballooning_grade = ballooning130,
  fibrosis_stage = fibrosis130,
  nas = nas130,
  stringsAsFactors = FALSE
)

counts130_names <- names(fread(
  p("data", "raw", "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz"),
  nrows = 0L, showProgress = FALSE
))[-1L]
stopifnot(setequal(counts130_names, meta130$title))
meta130 <- meta130[match(counts130_names, meta130$title), , drop = FALSE]
stopifnot(!anyNA(meta130$sample_id), identical(meta130$title, counts130_names))

p135 <- read_geo_metadata(p("data", "raw", "GSE135917_series_matrix.txt.gz"))
title135 <- tolower(as.character(p135$title))
study135 <- as.character(p135[["description.1"]])
group135 <- ifelse(grepl("normal subject", title135), "Control",
                   ifelse(grepl("before cpap", title135) |
                            (grepl("osa subject", title135) & !grepl("after cpap", title135)),
                          "OSA_baseline",
                          ifelse(grepl("after cpap", title135), "CPAP_post", "Other")))
patient135 <- ifelse(
  study135 == "STUDY GROUP 2",
  sub(".*CPAP\\s+", "", sub("[bc]$", "", as.character(p135$title), ignore.case = TRUE), ignore.case = TRUE),
  sub(".*Subject\\s+", "", as.character(p135$title), ignore.case = TRUE)
)
meta135 <- data.frame(
  sample_id = rownames(p135),
  title = as.character(p135$title),
  study_group = study135,
  group = group135,
  patient_id = patient135,
  age = as.numeric(p135[["age:ch1"]]),
  bmi = as.numeric(p135[["bmi:ch1"]]),
  sex = as.character(p135[["Sex:ch1"]]),
  platform = as.character(p135$platform_id),
  stringsAsFactors = FALSE
)

write_csv_atomic(meta126, p("data", "derived", "metadata_GSE126848.csv"))
write_csv_atomic(meta130, p("data", "derived", "metadata_GSE130970.csv"))
write_csv_atomic(meta135, p("data", "derived", "metadata_GSE135917.csv"))

contract <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  spec_version = read_spec()$spec_version,
  GSE126848 = list(
    n_samples = nrow(meta126),
    groups = as.list(table(meta126$group)),
    disease_labels = as.list(table(meta126$disease_label)),
    sex = as.list(table(meta126$sex)),
    count_rows_including_header = raw126_dims[1],
    count_columns_including_gene_id = ncol(raw126_head),
    sample_mapping = "exact match: count-file raw column ID to GEO description field"
  ),
  GSE130970 = list(
    n_samples = nrow(meta130),
    groups = as.list(table(meta130$group)),
    sex = as.list(table(meta130$sex)),
    control_rule = "steatosis=0 & fibrosis=0 & ballooning=0",
    counts_columns_match_titles = identical(meta130$title, counts130_names)
  ),
  GSE135917 = list(
    n_samples = nrow(meta135),
    groups = as.list(table(meta135$group)),
    study_groups = as.list(table(meta135$study_group)),
    discovery_groups = as.list(table(meta135$group[meta135$study_group == "STUDY GROUP 1"])),
    cpap_pairs = length(unique(meta135$patient_id[meta135$study_group == "STUDY GROUP 2"])),
    platform = unique(meta135$platform),
    covariates_complete = !anyNA(meta135[, c("age", "bmi", "sex")])
  )
)

stopifnot(
  nrow(meta126) == 57L,
  sum(meta126$group == "Control") == 26L,
  sum(meta126$group == "MASLD") == 31L,
  !anyNA(meta126$group),
  nrow(meta130) == 78L,
  sum(meta130$group == "Control") == 6L,
  sum(meta130$group == "MASLD") == 72L,
  !anyNA(meta130[, c("group", "sex", "age")]),
  sum(meta135$group == "Control") == 8L,
  sum(meta135$group == "OSA_baseline") == 34L
  ,sum(meta135$study_group == "STUDY GROUP 1") == 18L
  ,sum(meta135$study_group == "STUDY GROUP 2") == 48L
  ,sum(meta135$group == "OSA_baseline" & meta135$study_group == "STUDY GROUP 1") == 10L
  ,length(unique(meta135$patient_id[meta135$study_group == "STUDY GROUP 2"])) == 24L
  ,all(meta135$platform == "GPL6244")
  ,!anyNA(meta135[, c("age", "bmi", "sex")])
)

write_json_atomic(contract, p("data", "derived", "data_contract.json"))
save_session_info("00_audit_inputs")
cat("Input audit passed.\n")
