#!/bin/bash

# Define the target directory and the command executable
TARGET_DIR="/home/bopelage/Documents/heuristique v1.2/wcnf_files"
EXECUTABLE="/home/bopelage/Documents/OLLLP/bin/olllp"

# Check if the target directory actually exists before proceeding
if [ ! -d "$TARGET_DIR" ]; then
  echo "Error: Directory '$TARGET_DIR' not found."
  exit 1
fi

# Loop through all items in the target directory
for file in "$TARGET_DIR"/*; do
  # Check if the item is a regular file
  if [ -f "$file" ]; then
    echo "Running OLLLP on: $file"
    
    # Define the output filename (e.g., file.wcnf becomes file.wcnf.out)
    output_file="${file}.out"
    
    # Execute the command and save the output
    "$EXECUTABLE" -cpu-lim=1 "$file" > "$output_file" 2>&1
    
  fi
done

echo "All files processed! The outputs have been saved next to your original files."