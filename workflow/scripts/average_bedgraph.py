""" 
Creates average bedgraph file from list of bedgraph files
"""

import pandas as pd

bedgraphs = snakemake.input["bg"]
average_bedgraph = snakemake.output["average_bedgraph"]

# Read first bedgraph file
df = pd.read_csv(bedgraphs[0], 
                 sep="\t", 
                 header=None, 
                 skiprows=1, 
                 low_memory=False)

# Add value column (4th) of all other bedgraph files to df
for bedgraph in bedgraphs[1:]:
    df = pd.concat([df, pd.read_csv(bedgraph, 
                                    sep="\t", 
                                    header=None, 
                                    skiprows=1, 
                                    low_memory=False)[3]], axis=1)
    
# Calculate mean of 4th to last columns
df["mean"] = df.iloc[:, 3:].mean(axis=1)

# Only keep first 3 columns and mean column
df = df.iloc[:, :4]

# Write to file
df.to_csv(average_bedgraph, sep="\t", header=False, index=False)

