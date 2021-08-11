#!/usr/bin/env python3.8

# modules
import pandas as pd
import os
import sys


# input files
in_1 = sys.argv[1]
in_2 = sys.argv[2]
output_dir = sys.argv[3]
tax_level = sys.argv[4]

# create workspace
os.makedirs(output_dir, exist_ok=True)

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
df.to_csv(f"{output_dir}/RES.{tax_level}.combined_info.tsv", sep="\t", header=True, index=False)
