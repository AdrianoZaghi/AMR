folder="AMR_archive_paired_end"

# Iterate over each file in the folder
for file in "$folder"/*
do
    # Check if it's a file
    if [ -f "$file" ]; then
        # Compress the file into a .gz archive
        gzip "$file"
    fi
done
