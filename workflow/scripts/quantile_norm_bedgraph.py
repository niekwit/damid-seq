import logging
import pandas as pd

# Setup logging
log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(message)s",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
)


def quantile_normalize(df):
    """
    Perform quantile normalization on a DataFrame.

    Parameters:
    df (pandas.DataFrame): The input DataFrame to be quantile normalized.

    Returns:
    pandas.DataFrame: The quantile normalized DataFrame.
    """
    rank_mean = df.stack().groupby(df.rank(method="first").stack().astype(int)).mean()
    return df.rank(method="min").stack().astype(int).map(rank_mean).unstack()


bedgraphs = snakemake.input["bg"]
assert len(bedgraphs) > 0, "No bedgraphs found"
bgs = "\n".join(bedgraphs)
logging.info(f"Applying quantile normalisation to:\n{bgs}")

data = {}
for filename in bedgraphs:
    df = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "score"],
        engine="python",
        skiprows=1,
    )
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
    output_filename = filename.replace(
        "-norm.gatc.bedgraph", ".quantile-norm.gatc.bedgraph"
    )
    df = pd.DataFrame(probes, columns=["chrom", "start", "end"])
    df["score"] = normalized_scores[i]
    df["score"] = df["score"].round(3)
    df.to_csv(output_filename, sep="\t", header=False, index=False)
