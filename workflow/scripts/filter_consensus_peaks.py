import pandas as pd

# Load Snakemake variables
bed = snakemake.input["bed"]
peaks = snakemake.input["peaks"]
cs = snakemake.input["cs"]
ext_bed = snakemake.output[0]
k = int(snakemake.params["k"])
max_size = int(snakemake.params["max_size"])
sample = snakemake.wildcards["bg_sample"]

# Load chrom_sizes into dictionary
# This is needed to make sure that the extended regions 
# do not extend past the beginning/end of the chromosome
chrom_sizes = {}
with open(cs) as f:
    for line in f:
       (key, value) = line.split()
       chrom_sizes[key] = int(value)

# Load individual peaks into data frames
ind_peaks = []
for peak in peaks:
    df = pd.read_csv(peak, sep="\t", header=None, low_memory=False)
    ind_peaks.append(df)

# Counters for logging
extended_peaks = 0
skipped_peaks = 0
count = 1

with open(bed, "r") as bed_in, open(ext_bed, "w") as bed_out:
    while True:
        # For each line, reset extend_by value 
        # (may be set to 0 if region is too large)
        extend_by = int(snakemake.params["e"])
        
        # Read in one line at a time
        line = bed_in.readline()
        
        # Stop if end of file is reached
        if not line:
            break
        
        # Load line into variables
        chrom, start, end, num, *args = line.split()
        
        # Check how many regions are overlapping
        # Keep only regions that come from at least k regions
        if int(num) < k:
            skipped_peaks += 1
            continue # Skip to next line
        
        # Instead of only the overlapping region, get the region that contains 
        # the whole area of the overlapping regions
        starts = []
        ends = []
        for df in ind_peaks:
            df = df[(df[0] == chrom) & (df[1] <= int(end)) & (df[2] >= int(start))]
            if not df.empty:
                starts.append(df[1].min())
                ends.append(df[2].max())
        start = min(starts)
        end = max(ends)
                
        # Check if region is short enough to extend 
        # otherwise just keep the original region
        start = int(start)
        end = int(end)
        size = end - start
        if size >= max_size:
            extend_by = 0
                       
        # Extend the regions by x bp on each side (read from config.yaml)
        start = start - extend_by
        end = end + extend_by
        if extend_by != 0:
            extended_peaks += 1
        
        # Make sure that the extended regions do not 
        # extend past the begin/end of the chromosome
        if start < 0:
            start = 0
        if end > chrom_sizes[chrom]:
            end = chrom_sizes[chrom]
        
        # Generate name for extended region
        name = f"peak_{count}"
        count += 1
        
        # Write extended region to output file
        bed_out.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{num}\n")
        
with open(snakemake.log[0], "w") as log:
    log.write(f"Total number of peaks analysed: {count + skipped_peaks}\n")
    log.write(f"Skipped peaks: {skipped_peaks}\n")
    log.write(f"Overlapping peaks: {count}\n")
    log.write(f"Extended peaks: {extended_peaks}\n")
    
