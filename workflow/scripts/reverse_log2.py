import csv

# Get input/output files
bg_input = snakemake.input["bg"]
bg_output = snakemake.output["bg"]

with open(bg_input, "r") as infile, open(bg_output, "w") as outfile:
    # Read and write the file
    reader = csv.reader(infile, delimiter="\t")
    writer = csv.writer(outfile, delimiter="\t")

    # Skip the first row
    next(reader)

    # Iterate over the remaining rows to reverse the log2 transformation
    for row in reader:
        # Reverse the log2 transformation
        row[3] = 2 ** float(row[3])
        writer.writerow(row)
