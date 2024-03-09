import sys
import pandas as pd

def quantile_normalize(df):
    rank_mean = df.stack().groupby(df.rank(method="first").stack().astype(int)).mean()
    return df.rank(method="min").stack().astype(int).map(rank_mean).unstack()

if len(sys.argv) < 2:
    print("Please provide one or more .bedgraph files as command-line arguments.")
    sys.exit(1)

bedgraphs = sys.argv[1:]
bgs= "\n".join(bedgraphs)
print(f"Applying quantile normalisation to:\n{bgs}")

data = {}
for filename in bedgraphs:
    df = pd.read_csv(filename, 
                     sep="\t", 
                     header=None, 
                     names=["chrom", "start", "end", "score"], 
                     engine="python",
                     skiprows=1)
    df[["start", "end"]] = df[["start", "end"]].astype(int)
    for row in df.itertuples(index=False):
        key = (row.chrom, row.start)
        if key not in data:
            data[key] = [row.end, []]
        data[key][1].append(row.score)

probes = []
scores = []
for key, value in data.items():
    if len(value[1]) == len(bedgraphs):
        probes.append(list(key) + [value[0]])
        scores.append(value[1])

normalized_scores = quantile_normalize(pd.DataFrame(scores))

for i, filename in enumerate(bedgraphs):
    output_filename = filename.replace(".kde-norm.gatc.bedgraph", ".quantile-norm.gatc.bedgraph")
    df = pd.DataFrame(probes, columns=["chrom", "start", "end"])
    df["score"] = normalized_scores[i]
    df["score"] = df["score"].round(3)
    df.to_csv(output_filename, sep="\t", header=False, index=False)
