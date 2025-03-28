#!/usr/bin/env Rscript

# Read the GDS results file
gds_file <- "gds_result (3).txt"
gds_text <- readLines(gds_file)

# Initialize vectors to store data
entries <- list()

# Process each entry
current_entry <- list()
entry_number <- 0

# Debug: Print first few lines to check format
message("First few lines of input file:")
print(head(gds_text))

for (line in gds_text) {
  # Skip empty lines
  if (line == "") next
  
  # Check if this is a new entry (starts with a number and dot)
  if (grepl("^\\d+\\.", line)) {
    # If we have a previous entry, save it
    if (length(current_entry) > 0) {
      # Debug: Print current entry before saving
      message("\nProcessing entry:")
      print(current_entry)
      
      # Only save if we have all required fields
      if (!is.null(current_entry$accession) && 
          !is.null(current_entry$id) && 
          !is.null(current_entry$type) && 
          !is.null(current_entry$platform) && 
          !is.null(current_entry$ftp)) {
        entries[[length(entries) + 1]] <- current_entry
        message("Entry saved successfully")
      } else {
        message("Skipping incomplete entry: ", current_entry$accession)
        message("Missing fields:")
        if (is.null(current_entry$accession)) message("- accession")
        if (is.null(current_entry$id)) message("- id")
        if (is.null(current_entry$type)) message("- type")
        if (is.null(current_entry$platform)) message("- platform")
        if (is.null(current_entry$ftp)) message("- ftp")
      }
    }
    # Start new entry
    current_entry <- list()
    entry_number <- entry_number + 1
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

# Save the last entry if exists and complete
if (length(current_entry) > 0) {
  # Debug: Print last entry before saving
  message("\nProcessing last entry:")
  print(current_entry)
  
  if (!is.null(current_entry$accession) && 
      !is.null(current_entry$id) && 
      !is.null(current_entry$type) && 
      !is.null(current_entry$platform) && 
      !is.null(current_entry$ftp)) {
    entries[[length(entries) + 1]] <- current_entry
    message("Last entry saved successfully")
  } else {
    message("Skipping incomplete last entry: ", current_entry$accession)
    message("Missing fields:")
    if (is.null(current_entry$accession)) message("- accession")
    if (is.null(current_entry$id)) message("- id")
    if (is.null(current_entry$type)) message("- type")
    if (is.null(current_entry$platform)) message("- platform")
    if (is.null(current_entry$ftp)) message("- ftp")
  }
}

# Debug: Print number of entries collected
message("\nNumber of entries collected: ", length(entries))

# Create data frame from complete entries
if (length(entries) > 0) {
  gds_df <- do.call(rbind, lapply(entries, function(x) {
    data.frame(
      Accession = x$accession,
      ID = x$id,
      Type = x$type,
      Platform = x$platform,
      FTP_Download = x$ftp,
      stringsAsFactors = FALSE
    )
  }))
  
  # Save to CSV
  write.csv(gds_df, "gds_table_new.csv", row.names = FALSE)
  
  # Print summary
  message("\nCreated gds_table_new.csv with ", nrow(gds_df), " entries")
  message("Number of RNA-seq datasets: ", sum(grepl("RNA-seq|RNA sequencing", gds_df$Type, ignore.case = TRUE)))
  message("Number of microarray datasets: ", sum(grepl("array", gds_df$Type, ignore.case = TRUE)))
  
  # Debug: Print first few rows of the data frame
  message("\nFirst few rows of the data frame:")
  print(head(gds_df))
} else {
  message("No entries were collected. Check the input file format.")
} 