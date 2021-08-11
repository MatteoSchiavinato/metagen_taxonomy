#!/usr/bin/env python3.8

# modules
import pandas as pd
import os

# directories for raw data
raw_data_dir = "/gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/raw_data"
diversity_dir = "/gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/taxonomy/diversity"
output_dir = "/gpfs/data/fs71579/schmat90/CF/projects/meta/M00003_Toasted_soy_diversity/taxonomy/diversity/plots"

# create workspace
os.makedirs(output_dir, exist_ok=True)

# input files
in_1 = f"{raw_data_dir}/sample_description.tsv"
in_2 = f"{diversity_dir}/RES.G.counts.diversity.tsv"

# read information dataframe
df1 = pd.read_csv(in_1, sep="\t", comment="#", header=0, usecols=range(0,5), index_col=0)
df1.index = df1.index.rename("Sample")

# read diversity dataframe
df2 = pd.read_csv(in_2, sep="\t", comment="#", header=0)
df2.Sample = df2.Sample.str.replace("_.*", "", regex=True)
df2.index = df2.Sample
df2 = df2.iloc[:,1:]

# merge the two dataframes
df = df1.join(df2, how="left")
df = df.reset_index()

# write to output
df.to_csv(f"{output_dir}/RES.combined_info.tsv", sep="\t", header=True, index=False)
