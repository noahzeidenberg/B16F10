#!/usr/bin/env Rscript

# Read the GDS results file
gds_file <- "gds_result (3).txt"
gds_text <- readLines(gds_file)

# Initialize vectors to store data
accessions <- c()
ids <- c()
types <- c()
platforms <- c()
ftp_downloads <- c()

# Process each entry
current_entry <- list()
for (line in gds_text) {
  # Skip empty lines
  if (line == "") next
  
  # Check if this is a new entry (starts with a number and dot)
  if (grepl("^\\d+\\.", line)) {
    # If we have a previous entry, save it
    if (length(current_entry) > 0) {
      accessions <- c(accessions, current_entry$accession)
      ids <- c(ids, current_entry$id)
      types <- c(types, current_entry$type)
      platforms <- c(platforms, current_entry$platform)
      ftp_downloads <- c(ftp_downloads, current_entry$ftp)
    }
    # Start new entry
    current_entry <- list()
    next
  }
  
  # Extract information from each line
  if (grepl("^Type:", line)) {
    current_entry$type <- sub("^Type:\\s*", "", line)
  } else if (grepl("^Platform:", line)) {
    current_entry$platform <- sub("^Platform:\\s*", "", line)
  } else if (grepl("^FTP download:", line)) {
    current_entry$ftp <- sub("^FTP download:\\s*", "", line)
  } else if (grepl("^Series\\s+Accession:", line)) {
    current_entry$accession <- sub("^Series\\s+Accession:\\s*", "", line)
  } else if (grepl("^ID:", line)) {
    current_entry$id <- sub("^ID:\\s*", "", line)
  }
}

# Save the last entry if exists
if (length(current_entry) > 0) {
  accessions <- c(accessions, current_entry$accession)
  ids <- c(ids, current_entry$id)
  types <- c(types, current_entry$type)
  platforms <- c(platforms, current_entry$platform)
  ftp_downloads <- c(ftp_downloads, current_entry$ftp)
}

# Create data frame
gds_df <- data.frame(
  Accession = accessions,
  ID = ids,
  Type = types,
  Platform = platforms,
  FTP_Download = ftp_downloads
)

# Save to CSV
write.csv(gds_df, "gds_table_new.csv", row.names = FALSE)

# Print summary
message("Created gds_table_new.csv with ", nrow(gds_df), " entries")
message("Number of RNA-seq datasets: ", sum(grepl("RNA-seq|RNA sequencing", gds_df$Type, ignore.case = TRUE))) 