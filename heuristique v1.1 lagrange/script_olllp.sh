#!/bin/bash

# Define the target directory and the command executable
TARGET_DIR=$1
#TARGET_DIR="./wcnf_files"
EXECUTABLE="/home/bopelage/bin/olllp"

echo "Problem,Time,Cost" > "$2"

# Check if the target directory actually exists before proceeding
if [ ! -d "$TARGET_DIR" ]; then
  echo "Error: Directory '$TARGET_DIR' not found."
  exit 1
fi

# Loop through all items in the target directory
for file in "$TARGET_DIR"/*.wcnf; do
  # Check if the item is a regular file
  if [ -f "$file" ]; then
    echo "Running OLLLP on: $file"
    
    # Define the output filename (e.g., file.wcnf becomes file.wcnf.out)
    output_file="${file%.*}_${3}.out"
    
    # Execute the command and save the output
    "$EXECUTABLE" -cpu-lim=${3} "$file" > "$output_file" 2>&1
    awk 'BEGIN { time=0; cost=0 } \
         /INTERRUPTED/ { time='$3' } \
         /c CPU:/ { time=0+$3 } \
         /o / { cost=$2 } \
         END { print "'$file'," time "," cost }' "$output_file" >> "$2"
  fi
done

echo "All files processed! The outputs have been saved next to your original files."