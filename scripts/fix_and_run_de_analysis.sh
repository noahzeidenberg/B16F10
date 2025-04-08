#!/bin/bash
# Script to fix the mapping and design matrix issues and run the differential expression analysis

# Set the GSE ID
GSE_ID="GSE287957"

# Set the base directory
BASE_DIR="$HOME/scratch/B16F10"

# Create logs directory if it doesn't exist
mkdir -p "$BASE_DIR/logs/fix_de_analysis"

# Run the mapping file creation script
echo "Creating mapping file..."
Rscript "$BASE_DIR/scripts/create_mapping_file.R" > "$BASE_DIR/logs/fix_de_analysis/create_mapping_file.log" 2>&1

# Check if the script ran successfully
if [ $? -ne 0 ]; then
  echo "Error creating mapping file. Check the log file for details."
  exit 1
fi

# Run the design matrix update script
echo "Updating design matrix..."
Rscript "$BASE_DIR/scripts/update_design_matrix.R" > "$BASE_DIR/logs/fix_de_analysis/update_design_matrix.log" 2>&1

# Check if the script ran successfully
if [ $? -ne 0 ]; then
  echo "Error updating design matrix. Check the log file for details."
  exit 1
fi

# Run the preprocessing script
echo "Preprocessing normalized counts..."
Rscript "$BASE_DIR/scripts/preprocess_counts.R" > "$BASE_DIR/logs/fix_de_analysis/preprocess_counts.log" 2>&1

# Check if the script ran successfully
if [ $? -ne 0 ]; then
  echo "Error preprocessing normalized counts. Check the log file for details."
  exit 1
fi

# Run the fix design matrix script
echo "Fixing design matrix..."
Rscript "$BASE_DIR/scripts/fix_design_matrix.R" "$GSE_ID" > "$BASE_DIR/logs/fix_de_analysis/fix_design_matrix.log" 2>&1

# Check if the script ran successfully
if [ $? -ne 0 ]; then
  echo "Error fixing design matrix. Check the log file for details."
  exit 1
fi

# Run the differential expression analysis
echo "Running differential expression analysis..."
Rscript "$BASE_DIR/scripts/differential_expression.R" "$GSE_ID" > "$BASE_DIR/logs/fix_de_analysis/differential_expression.log" 2>&1

# Check if the script ran successfully
if [ $? -ne 0 ]; then
  echo "Error running differential expression analysis. Trying modified approach..."
  
  # Run the modified differential expression analysis
  echo "Running modified differential expression analysis..."
  Rscript "$BASE_DIR/scripts/differential_expression_modified.R" "$GSE_ID" > "$BASE_DIR/logs/fix_de_analysis/differential_expression_modified.log" 2>&1
  
  # Check if the script ran successfully
  if [ $? -ne 0 ]; then
    echo "Error running modified differential expression analysis. Trying fixed approach..."
    
    # Run the fixed differential expression analysis
    echo "Running fixed differential expression analysis..."
    Rscript "$BASE_DIR/scripts/differential_expression_fixed.R" "$GSE_ID" > "$BASE_DIR/logs/fix_de_analysis/differential_expression_fixed.log" 2>&1
    
    # Check if the script ran successfully
    if [ $? -ne 0 ]; then
      echo "Error running fixed differential expression analysis. Check the log file for details."
      exit 1
    fi
    
    echo "Fixed differential expression analysis completed successfully."
  else
    echo "Modified differential expression analysis completed successfully."
  fi
else
  echo "Differential expression analysis completed successfully."
fi

echo "All steps completed successfully." 