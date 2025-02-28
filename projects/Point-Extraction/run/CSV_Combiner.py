import glob

# List of input file names
input_files = glob.glob("*.csv")  # Use the appropriate pattern to match your input files

# Output file name
output_file = "combined_output.csv"

# Write header to the output file
with open(output_file, "w") as output_csv:
    output_csv.write("GP_x0,GP_y0,GP_z0\n")

    # Loop through each input file
    for file_name in input_files:
        with open(file_name, "r") as input_csv:
            # Skip the header line in each input file
            next(input_csv)
            # Copy the rest of the lines to the output file
            for line in input_csv:
                output_csv.write(line)

print("Files combined successfully. Output written to:", output_file)
