#!/usr/bin/env python3

# modules
import pandas as pd
import numpy as np
import sys
import os
import argparse as ap
from time import asctime as at

# arguments
p = ap.ArgumentParser()
p.add_argument("--input-dir", help="Input directory containing bracken report files (use only if not using --input-files)")
p.add_argument("--extension", help="Extension of bracken report files contained in --input-dir (default: \"report\")", type=str, default="report")
p.add_argument("--input-files", help="Input files to be combined (use only if not using --input-dir and --extension)", nargs="+")
p.add_argument("--classif-level", help="Either one of S, G, F, C, O, P, D", default="S", type=str)
p.add_argument("--max-species", help="Maximum number of species to include in final table, determined by highest average frequency", type=int, default=10)
p.add_argument("--output-dir", help="Directory where to put output files", default=".")
args = p.parse_args()

# output workspace
os.makedirs(args.output_dir, exist_ok=True)

# conditions
if args.max_species > 20:
	sys.stderr.write("\n### ERROR: don't select more than 20 max species with --max-species, this will complicate the plotting and is not implemented.\n\n")

# read input files
sys.stderr.write("[{0}] Reading input file paths\n".format(at()))
if args.input_dir and args.input_files:
	sys.stderr.write("\n### ERROR: can't use both --input-dir and --input-files\n\n")
	sys.exit()

if args.input_dir:
	Input_files = [ args.input_dir + "/" + i for i in os.listdir(args.input_dir) if i.split(".")[-1] == args.extension ]
elif args.input_files:
	Input_files = args.input_files

# create table
sys.stderr.write("[{0}] Creating empty dataframes\n".format(at()))
df_frac = pd.DataFrame()
df_count = pd.DataFrame()

# iterate over all files
sys.stderr.write("[{0}] Filling empty dataframes with actual data from input\n".format(at()))
for FILENAME in Input_files:

	# read dataframe
	tmp = pd.read_csv(FILENAME, sep="\t", index_col=5, header=0)

	# extract sample name
	sample = FILENAME.split("/")[-1].split(".")[0]

	# rename columns and index
	tmp.columns = ["Fraction", "Counts", "Added", "Taxrank", "Taxid"]
	tmp.index = tmp.index.rename("Taxa")

	# select only lines corresponding to the specific classification level
	mask = tmp["Taxrank"] == args.classif_level
	tmp = tmp.loc[mask , :]

	# remove whitespaces
	tmp.index = tmp.index.str.strip(" ")

	# subdivide into fraction and counts dataframes
	tmp_frac = tmp.loc[: , ["Fraction"]]
	tmp_count = tmp.loc[: , ["Counts"]]

	# rename column with sample name
	tmp_frac.columns = [sample]
	tmp_count.columns = [sample]

	# join with main dataframe
	df_frac = df_frac.join(tmp_frac, how="outer")
	df_count = df_count.join(tmp_count, how="outer")

# sorting dataframe columns
sys.stderr.write("[{0}] Sorting dataframe\n".format(at()))
df_frac.columns = sorted(df_frac.columns.to_list())
df_count.columns = sorted(df_count.columns.to_list())

# remove NaN
sys.stderr.write("[{0}] Converting NaN to 0\n".format(at()))
df_frac = df_frac.fillna(0)
df_count = df_count.fillna(0)

# sum duplicates
df_frac = df_frac.groupby(df_frac.index).sum()
df_count = df_count.groupby(df_count.index).sum()

# create subset table with top species
Top_entries = df_frac.mean(axis=1).sort_values(ascending=False).head(args.max_species).index.to_list()
Nontop_entries = [entry for entry in df_frac.index.to_list() if entry not in Top_entries]
df_frac_nontop = df_frac.loc[Nontop_entries , :]
df_top = df_frac.loc[Top_entries , :]
df_nontop = df_frac_nontop.sum(axis=0)


df_nontop = df_nontop.rename("Other")
df_top = df_top.append(df_nontop)

# rescale values to 100% (to account for fluctuations)
# kraken and bracken not always produce exact 100% totals
# this applies only to the df_frac database of course
sys.stderr.write("[{0}] Rescaling relative abundances to 100%\n".format(at()))
df_frac = df_frac / df_frac.sum() * 100
df_top = df_top / df_top.sum() * 100

# reset index
sys.stderr.write("[{0}] Preparing for output\n".format(at()))
df_frac = df_frac.reset_index()
df_top = df_top.reset_index()
df_count = df_count.reset_index()

# write to output
sys.stderr.write("[{0}] Writing to output\n".format(at()))
outfile_frac = args.output_dir + "/" + ".".join(["RES", args.classif_level, "rel_abundance", "tsv"])
outfile_top = args.output_dir + "/" + ".".join(["RES", args.classif_level, "rel_abundance", "top_" + str(args.max_species), "tsv"])
outfile_count = args.output_dir + "/" + ".".join(["RES", args.classif_level, "counts", "tsv"])
df_frac.to_csv(outfile_frac, sep="\t", header=True, index=False)
df_top.to_csv(outfile_top, sep="\t", header=True, index=False)
df_count.to_csv(outfile_count, sep="\t", header=True, index=False)
