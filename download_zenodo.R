cat("All extractions complete!\n")\ncat("Files saved to:", normalizePath(output_dir), "\n")7"
output_dir <- "data/raw"
zenodo_api <- "https://zenodo.org/api/records"

# Create output directory if it doesn't exist
dir_create(output_dir, recurse = TRUE)

# Extract record ID from DOI
# Zenodo DOI format: 10.5281/zenodo.RECORD_ID
record_id <- sub("^.*\\.", "", doi)

cat("Downloading files from Zenodo record:", record_id, "\n")

# Fetch metadata from Zenodo API
api_url <- paste0(zenodo_api, "/", record_id)
response <- GET(api_url)

if (status_code(response) != 200) {
  stop("Failed to retrieve record from Zenodo. Status code: ", status_code(response))
}

# Parse JSON response
record_data <- fromJSON(content(response, as = "text"))

# Extract files
files <- record_data$files$entries
n_files <- length(files)

cat("Found", n_files, "file(s) to download\n\n")

# Download each file
for (i in seq_along(files)) {
  file_name <- files[[i]]$key
  file_url <- files[[i]]$links$self
  file_size <- files[[i]]$size
  
  output_path <- file.path(output_dir, file_name)
  
  cat("[", i, "/", n_files, "] Downloading:", file_name, 
      "(", round(file_size / 1e6, 2), "MB )\n", sep = "")
  
  # Download file with progress
  tryCatch(
    {
      download.file(
        url = file_url,
        destfile = output_path,
        mode = "wb",
        quiet = FALSE
      )
      cat("✓ Successfully downloaded to:", output_path, "\n\n")
    },
    error = function(e) {
      cat("✗ Error downloading:", file_name, "\n")
      cat("  Error:", e$message, "\n\n")
    }
  )
}

cat("Download complete!\n\n")

# Unzip NTA files into organized subdirectories
cat("Processing NTA zip files...\n")

nta_base_dir <- file.path(output_dir, "NTA")
nta_dg_dir <- file.path(nta_base_dir, "DG")
nta_no_dg_dir <- file.path(nta_base_dir, "No_DG")

# Create NTA subdirectories
dir_create(nta_dg_dir, recurse = TRUE)
dir_create(nta_no_dg_dir, recurse = TRUE)

# Process NTA_Density_gradient.zip
dg_zip <- file.path(output_dir, "NTA_Density_gradient.zip")
if (file_exists(dg_zip)) {
  cat("Extracting NTA_Density_gradient.zip to", nta_dg_dir, "...\n")
  tryCatch(
    {
      unzip(dg_zip, exdir = nta_dg_dir)
      cat("✓ Successfully extracted NTA_Density_gradient.zip\n\n")
    },
    error = function(e) {
      cat("✗ Error extracting NTA_Density_gradient.zip:", e$message, "\n\n")
    }
  )
} else {
  cat("⚠ NTA_Density_gradient.zip not found in", output_dir, "\n\n")
}

# Process NTA_NoDensityGradient.zip
no_dg_zip <- file.path(output_dir, "NTA_NoDensityGradient.zip")
if (file_exists(no_dg_zip)) {
  cat("Extracting NTA_NoDensityGradient.zip to", nta_no_dg_dir, "...\n")
  tryCatch(
    {
      unzip(no_dg_zip, exdir = nta_no_dg_dir)
      cat("✓ Successfully extracted NTA_NoDensityGradient.zip\n\n")
    },
    error = function(e) {
      cat("✗ Error extracting NTA_NoDensityGradient.zip:", e$message, "\n\n")
    }
  )
} else {
  cat("⚠ NTA_NoDensityGradient.zip not found in", output_dir, "\n\n")
}

cat("Processing complete!\n")
cat("Files saved to:", normalizePath(output_dir), "\n")
cat("NTA files organized in:", normalizePath(nta_base_dir), "\n")")



cat("Processing complete!\n")
cat("Files saved to:", normalizePath(output_dir), "\n")
cat("NTA files organized in:", normalizePath(nta_base_dir), "\n")

# Process FCS_files.zip
cat("\nProcessing FCS files...\n")

fcm_dir <- file.path(output_dir, "FCM")
dir_create(fcm_dir, recurse = TRUE)

fcs_zip <- file.path(output_dir, "FCS_files.zip")
if (file_exists(fcs_zip)) {
  cat("Extracting FCS_files.zip to", fcm_dir, "...\n")
  tryCatch(
    {
      unzip(fcs_zip, exdir = fcm_dir)
      cat("✓ Successfully extracted FCS_files.zip\n\n")
    },
    error = function(e) {
      cat("✗ Error extracting FCS_files.zip:", e$message, "\n\n")
    }
  )
} else {
  cat("⚠ FCS_files.zip not found in", output_dir, "\n\n")
}

# Process PD_output.zip
cat("Processing PD_output files...\n")

pd_zip <- file.path(output_dir, "PD_output.zip")
if (file_exists(pd_zip)) {
  cat("Extracting PD_output.zip to", output_dir, "...\n")
  tryCatch(
    {
      unzip(pd_zip, exdir = output_dir)
      cat("✓ Successfully extracted PD_output.zip\n\n")
    },
    error = function(e) {
      cat("✗ Error extracting PD_output.zip:", e$message, "\n\n")
    }
  )
} else {
  cat("⚠ PD_output.zip not found in", output_dir, "\n\n")
}

# Process detected_proteins_annotation.zip
cat("Processing detected proteins annotation files...\n")

proteins_zip <- file.path(output_dir, "detected_proteins_annotation.zip")
if (file_exists(proteins_zip)) {
  cat("Extracting detected_proteins_annotation.zip to", output_dir, "...\n")
  tryCatch(
    {
      unzip(proteins_zip, exdir = output_dir)
      cat("✓ Successfully extracted detected_proteins_annotation.zip\n\n")
    },
    error = function(e) {
      cat("✗ Error extracting detected_proteins_annotation.zip:", e$message, "\n\n")
    }
  )
} else {
  cat("⚠ detected_proteins_annotation.zip not found in", output_dir, "\n\n")
}

cat("All processing complete!\n")
cat("FCM files organized in:", normalizePath(fcm_dir), "\n")

cat("All extractions complete!\n")
cat("Files saved to:", normalizePath(output_dir), "\n")")