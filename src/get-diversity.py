#!/usr/bin/env python3

# modules
import pandas as pd
import numpy as np
import sys
from plotnine import *
from skbio import diversity as biodiv
from skbio.stats.distance import DistanceMatrix, permanova
from scipy.spatial.distance import pdist, squareform
from skbio.stats.ordination import pcoa, OrdinationResults
from sklearn import manifold
from sklearn import preprocessing
import de_toolkit as detk
import re
from itertools import combinations
import argparse as ap

# arguments
p = ap.ArgumentParser()
p.add_argument("--input-file", help="Matrix of counts with samples as columns and IDs as rows")
args = p.parse_args()

############
### TAXA ###
############


# regenerate dictionaries
Timepoints = ["D3", "D14", "D21", "D35"]
Treatments = ["CON", "PFA", "AGP", "AGP+PFA", "AB", "AB+PFA"]
Replicates = ["1", "2", "3", "4"]


###########################
### ALPHA-DIVERSITY OF TAXA

input_file = "/binfl/lv70694/schmat/chicken/bracken/RESULTS.genus.counts.txt"
df = pd.read_csv(input_file, sep="\t")
df = df.set_index("Taxa")

# transposing the matrix
df = df.T
df.index = df.index.str.replace("^F", "D")

# pcoa
# 'eigh' = exact eigenvectors and eigenvalues for all dimensions
names = df.index
x = squareform(pdist(df, metric="braycurtis"))
x = DistanceMatrix(x, ids=names)
x = pcoa(x, method='eigh')

comp = x.samples.loc[:,["PC1","PC2","PC3"]]
comp.index = names
comp["Treatment"] = comp.index.str.split("_").str[1]
comp["Treatment"] = pd.Categorical(comp.Treatment, categories=Treatments)
comp["Timepoint"] = comp.index.str[0:3].str.replace("_", "").str.replace("^F", "D")
comp["Timepoint"] = pd.Categorical(comp.Timepoint, categories=Timepoints)

comp_pcoa = comp.copy()
comp_pcoa["Code"] = comp_pcoa["Timepoint"].astype("str") + "_" + comp_pcoa["Treatment"].astype("str")
comp_pcoa.Code = comp_pcoa.Code.str.replace("^D3_.*", "D3")
comp_pcoa = comp_pcoa.loc[:,["PC1", "PC2", "PC3", "Code"]]

subCategories = ["D14_AGP+PFA", "D14_AB+PFA", "D21_AGP+PFA", "D21_AB+PFA", "D35_AGP+PFA", "D35_AB+PFA"]
mask = [i in subCategories for i in comp_pcoa.Code.tolist()]
comp_pcoa_sub = comp_pcoa.loc[mask , :]
comp_pcoa_sub["Code"] = pd.Categorical(comp_pcoa_sub.Code, categories=subCategories)


# plot

### pcoa with ALL samples ###

comp_pcoa.Code = pd.Categorical(comp_pcoa.Code, categories=myColors.keys())

P1 = (
ggplot(data=comp_pcoa)
+ geom_point(mapping=aes(x="PC1", y="PC2", color="Code", fill="Code"), size=3, show_legend=None)
+ ggtitle("")
+ xlab("PC1")
+ ylab("PC2")
+ scale_fill_manual(name="", values=myColors)
+ scale_color_manual(name="", values=myColors, guide=False)
+ scale_size_continuous(guide=False)
+ theme(plot_title=element_text(family="sans", face="bold",
					colour='#261A1A', size=12),
	axis_title_x=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_title_y=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_text_x=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	axis_text_y=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	legend_text=element_text(family="sans", face="plain", colour='#000000', size=10),
	legend_title=element_text(family="sans", face="bold", colour='#000000', size=10),
	panel_background=element_rect(fill="white", color="white"),
	panel_grid=element_line(color="#E6E6E6"),
	aspect_ratio=0.3,
	panel_spacing=0.5)
)

P1.save(filename="/binfl/lv70694/schmat/chicken/diversity_taxa/RESULTS.PCoA.taxa.1_vs_2.png",
format="png", height=15, width=25, dpi=300, units="cm")

P2 = (
ggplot(data=comp_pcoa)
+ geom_point(mapping=aes(x="PC1", y="PC3", color="Code", fill="Code"), size=3, show_legend=None)
+ ggtitle("")
+ xlab("PC1")
+ ylab("PC3")
+ scale_fill_manual(name="", values=myColors)
+ scale_color_manual(name="", values=myColors, guide=False)
+ scale_size_continuous(guide=False)
+ theme(plot_title=element_text(family="sans", face="bold",
					colour='#261A1A', size=12),
	axis_title_x=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_title_y=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_text_x=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	axis_text_y=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	legend_text=element_text(family="sans", face="plain", colour='#000000', size=10),
	legend_title=element_text(family="sans", face="bold", colour='#000000', size=10),
	panel_background=element_rect(fill="white", color="white"),
	panel_grid=element_line(color="#E6E6E6"),
	aspect_ratio=0.3,
	panel_spacing=0.5)
)

P2.save(filename="/binfl/lv70694/schmat/chicken/diversity_taxa/RESULTS.PCoA.taxa.1_vs_3.png",
format="png", height=15, width=25, dpi=300, units="cm")

P3 = (
ggplot(data=comp_pcoa)
+ geom_point(mapping=aes(x="PC2", y="PC3", color="Code", fill="Code"), size=3, show_legend=None)
+ ggtitle("")
+ xlab("PC2")
+ ylab("PC3")
+ scale_fill_manual(name="", values=myColors)
+ scale_color_manual(name="", values=myColors, guide=False)
+ scale_size_continuous(guide=False)
+ theme(plot_title=element_text(family="sans", face="bold",
					colour='#261A1A', size=12),
	axis_title_x=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_title_y=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_text_x=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	axis_text_y=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	legend_text=element_text(family="sans", face="plain", colour='#000000', size=10),
	legend_title=element_text(family="sans", face="bold", colour='#000000', size=10),
	panel_background=element_rect(fill="white", color="white"),
	panel_grid=element_line(color="#E6E6E6"),
	aspect_ratio=0.3,
	panel_spacing=0.5)
)

P3.save(filename="/binfl/lv70694/schmat/chicken/diversity_taxa/RESULTS.PCoA.taxa.2_vs_3.png",
format="png", height=15, width=25, dpi=300, units="cm")


# calculating alpha diversity
# Species richness is the number of different species in a sample.
# Shannon index measures how evenly the microbes are distributed in a sample.
richness = biodiv.alpha_diversity('observed_otus', df, df.index, validate=True)
diversity = biodiv.alpha_diversity('shannon', df, df.index, validate=True)

# get dataframe for plotting
timepoint = df.index.str.split("_").str[0]
timepoint = pd.Categorical(timepoint, categories=timepoint.unique())
treatment = df.index.str.split("_").str[1]
treatment = pd.Categorical(treatment, categories=treatment.unique())

x = pd.DataFrame({"Timepoint":timepoint, "Treatment":treatment, "Richness":richness})
y = pd.DataFrame({"Timepoint":timepoint, "Treatment":treatment, "Diversity":diversity})

# replace D3 treatments with "D3" (they're all the same)
x["Treatment"] = x["Treatment"].astype("object")
x.loc[x["Timepoint"]=="D3" , ["Treatment"]] = "D3"

y["Treatment"] = y["Treatment"].astype("object")
y.loc[x["Timepoint"]=="D3" , ["Treatment"]] = "D3"

# set categorical
x["Treatment"] = pd.Categorical(x["Treatment"], categories=x["Treatment"].unique())
x["Timepoint"] = pd.Categorical(x["Timepoint"], categories=x["Timepoint"].unique())

y["Treatment"] = pd.Categorical(y["Treatment"], categories=y["Treatment"].unique())
y["Timepoint"] = pd.Categorical(y["Timepoint"], categories=y["Timepoint"].unique())

# rename D3
newIndex = x.index.tolist()
newIndex[0:6] = ["D3_CON_R1", "D3_CON_R2", "D3_CON_R3", "D3_CON_R4", "D3_CON_R5", "D3_CON_R6"]
x.index = newIndex
y.index = newIndex

# write to output the data tables
x.to_csv("/binfl/lv70694/schmat/chicken/diversity_taxa/RESULTS.richness.all.txt", sep="\t", header=True, index=True)
y.to_csv("/binfl/lv70694/schmat/chicken/diversity_taxa/RESULTS.diversity.all.txt", sep="\t", header=True, index=True)


###################################
### BETA-DIVERSITY OF ARG GENES ###
###################################

myColors={
"D3":"#242424",
"D14_CON":"#C6C6C6",
"D14_PFA":"#FFB788",
"D14_AGP":"#52B788",
"D14_AGP+PFA":"#EEAEEE",
"D14_AB":"#8CD9FF",
"D14_AB+PFA":"#EDE2AD",
"D21_CON":"#979797",
"D21_PFA":"#FF805F",
"D21_AGP":"#38805F",
"D21_AGP+PFA":"#FF71F1",
"D21_AB":"#80A1F5",
"D21_AB+PFA":"#FFD970",
"D35_CON":"#575757",
"D35_PFA":"#FF5949",
"D35_AGP":"#365949",
"D35_AGP+PFA":"#943B8B",
"D35_AB":"#3B4F95",
"D35_AB+PFA":"#CDA153",
}

TimeColors = {"D3":"#808080", "D14":"#1FE5F7", "D21":"#D5894A", "D35":"#3BAA3B"}
TreatColors = {"D3":"#808080", "CON": "#FF934F", "PFA":"#EF6461", "AGP":"#E4B363",
				"AGP+PFA":"#AD7A9A", "AB":"#6A8D73", "AB+PFA":"#35B1E6"}


# Bray–Curtis dissimilarity
# - based on abundance or read count data
# - differences in microbial abundances between two samples (e.g., at species level)
#     values are from 0 to 1
#     0 means both samples share the same species at exact the same abundances
#     1 means both samples have complete different species abundances

# Bray-Curtis dissimilarity calculation, which has a number of ideal properties:
#
#     It is invariant to changes in units.
#     It is unaffected by additions/removals of species that are not present in two communities,
#     It is unaffected by the addition of a new community,
#     It can recognize differences in total abundances when relative abundances are the same,

# remove D3
df = df[~df.index.str.contains("D3_")]

grouping = [ i.split("_")[0:2] for i in df.index.tolist() ]
grouping = [ "_".join(lst) + "_" + "R" for lst in grouping ]
# grouping = [ re.sub("D3_.*", "D3_CON_R", i) for i in grouping ]
pairs = list(set(grouping))
Combinations = combinations(pairs, 2)
# Combinations = [ i for i in Combinations if (i[0][0:3] == i[1][0:3]) or (i[0]=="D3_CON_R") or (i[1] == "D3_CON_R") ]
Combinations = [ i for i in Combinations if (i[0][0:3] == i[1][0:3]) ]

# permanova

# for pair in Combinations
# get lines from df that match the two pairs
# if D3_CON: special condition
# assign grouping
# perform permanova
# append p.value and comparison to a list
# write list to file

PNOVA_OUT = open("/binfl/lv70694/schmat/chicken/diversity_taxa/RESULTS.dissimilarity.permanova.txt", "w")
PNOVA_OUT.write("\t".join(["Group_1", "Group_2", "P_value", "<0.05"]) + "\n")
for pair in Combinations:

	couple = list(pair)
	#
	# if couple[0] == "D3_CON_R":
	# 	couple[0] = "D3_"
	# if couple[1] == "D3_CON_R":
	# 	couple[1] = "D3_"

	mask = [i for i in df.index.tolist() if (couple[0] in i) or (couple[1] in i)]
	tmp = df.loc[mask , :]
	tmp_grouping = [ i.split("_")[0:2] for i in tmp.index.tolist() ]
	tmp_grouping = [ "_".join(lst) for lst in tmp_grouping ]
	# tmp_grouping = [ re.sub("D3_.*", "D3", i) for i in tmp_grouping ]
	diss = biodiv.beta_diversity('braycurtis', tmp, tmp.index, validate=True)
	p_nova = permanova(diss, tmp_grouping, permutations=999)

	couple[0] = couple[0].replace("_R", "")
	couple[1] = couple[1].replace("_R", "")
	# if couple[0] == "D3_":
	# 	couple[0] = "D3"
	# if couple[1] == "D3_":
	# 	couple[1] = "D3"

	if p_nova["p-value"] < 0.05:
		significance = "PASS"
	else:
		significance = "FAIL"
	PNOVA_OUT.write("{0}\t{1}\t{2}\t{3}\n".format(couple[0], couple[1], p_nova["p-value"], significance))

PNOVA_OUT.close()

# distance matrix
dissimilarity = biodiv.beta_diversity('braycurtis', df, df.index, validate=True)
diss_df = pd.DataFrame(dissimilarity.data, index=dissimilarity.ids, columns=dissimilarity.ids)
embedding = manifold.MDS(n_components=2, dissimilarity="precomputed", metric=False)
diss_NMDS = embedding.fit_transform(diss_df)
diss_NMDS_df = pd.DataFrame(diss_NMDS, index=diss_df.index, columns=["X", "Y"])

# timepoint = df.index.str[0:3].str.replace("_", "").str.replace("^F", "D")
timepoint = df.index.str[0:3].str.replace("_", "")
timepoint = pd.Categorical(timepoint, categories=timepoint.unique())
treatment = df.index.str.split("_").str[1]
treatment = pd.Categorical(treatment, categories=treatment.unique())
x = pd.DataFrame({"Timepoint":timepoint, "Treatment":treatment, "X":diss_NMDS_df.X, "Y":diss_NMDS_df.Y})

# replace D3 treatments with "D3"
x["Treatment"] = x["Treatment"].astype("object")
x["Timepoint"] = x["Timepoint"].astype("object")
# x.loc[x.Timepoint == "D3", "Treatment"] = "D3"
x_sub = x.copy()
x_time = x.copy()
x_treat = x.copy()

objs = [x["Timepoint"] + "_" + x["Treatment"], x["Timepoint"], x["X"], x["Y"]]
x = pd.concat(objs, axis=1)
x.columns = ["Condition", "Timepoint", "X", "Y"]
# x.loc[x.Condition == "D3_D3" , ["Condition"]] = "D3"
x["Condition"] = x["Condition"].astype("object")

# extract meaningful comparison
# x_sub = x_sub.loc[(x_sub.Treatment == "D3") | (x_sub.Treatment == "AGP+PFA") | (x_sub.Treatment == "AB+PFA"), :]
x_sub = x_sub.loc[(x_sub.Treatment == "AGP+PFA") | (x_sub.Treatment == "AB+PFA"), :]
objs = [x_sub["Timepoint"] + "_" + x_sub["Treatment"], x_sub["Timepoint"], x_sub["X"], x_sub["Y"]]
x_sub = pd.concat(objs, axis=1)
x_sub.columns = ["Condition", "Timepoint", "X", "Y"]
# x_sub.loc[x_sub.Condition == "D3_D3" , ["Condition"]] = "D3"
x_sub["Condition"] = x_sub["Condition"].astype("object")

# nmds with time only (no treat)
# and with treat only (no time)
x_time = x_time.loc[:,["Timepoint", "X", "Y"]]
x_treat = x_treat.loc[:,["Treatment", "X", "Y"]]

# write to output
x.to_csv("/binfl/lv70694/schmat/chicken/diversity_taxa/RESULTS.dissimilarity.txt", sep="\t", header=True, index=True)
x_sub.to_csv("/binfl/lv70694/schmat/chicken/diversity_taxa/RESULTS.dissimilarity.D_vs_F.txt", sep="\t", header=True, index=True)

# set categories

# Categories=["D3",
# 			"D14_CON", "D21_CON", "D35_CON",
# 			"D14_PFA", "D21_PFA", "D35_PFA",
# 			"D14_AGP", "D21_AGP", "D35_AGP",
# 			"D14_AGP+PFA", "D21_AGP+PFA", "D35_AGP+PFA",
# 			"D14_AB", "D21_AB", "D35_AB",
# 			"D14_AB+PFA", "D21_AB+PFA", "D35_AB+PFA", ]

Categories=["D14_CON", "D21_CON", "D35_CON",
			"D14_PFA", "D21_PFA", "D35_PFA",
			"D14_AGP", "D21_AGP", "D35_AGP",
			"D14_AGP+PFA", "D21_AGP+PFA", "D35_AGP+PFA",
			"D14_AB", "D21_AB", "D35_AB",
			"D14_AB+PFA", "D21_AB+PFA", "D35_AB+PFA", ]

x["Condition"] = pd.Categorical(x.Condition, categories=Categories)

# Categories=["D3",
# 			"D14_AGP+PFA", "D21_AGP+PFA", "D35_AGP+PFA",
# 			"D14_AB+PFA", "D21_AB+PFA", "D35_AB+PFA", ]

Categories=["D14_AGP+PFA", "D21_AGP+PFA", "D35_AGP+PFA",
			"D14_AB+PFA", "D21_AB+PFA", "D35_AB+PFA", ]

x_sub["Condition"] = pd.Categorical(x_sub.Condition, categories=Categories)

# plot all
P8 = (
ggplot(data=x)
+ geom_point(mapping=aes(x="X", y="Y", color="Condition", fill="Condition"), size=3, show_legend=None)
+ stat_ellipse(	mapping=aes(x="X", y="Y", color="Condition"), size=0.6, show_legend=None,
				type="t", level=0.95, segments=51)
+ ggtitle("")
+ xlab("NMDS1")
+ ylab("NMDS2")
+ scale_fill_manual(name="", values=myColors)
+ scale_color_manual(name="", values=myColors, guide=False)
+ scale_size_continuous(guide=False)
+ theme(plot_title=element_text(family="sans", face="bold",
					colour='#261A1A', size=12),
	axis_title_x=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_title_y=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_text_x=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	axis_text_y=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	legend_text=element_text(family="sans", face="plain", colour='#000000', size=10),
	legend_title=element_text(family="sans", face="bold", colour='#000000', size=10),
	panel_background=element_rect(fill="white", color="white"),
	panel_grid=element_line(color="#E6E6E6"),
	aspect_ratio=0.3,
	panel_spacing=0.5)
)

P8.save(filename="/binfl/lv70694/schmat/chicken/diversity_taxa/RESULTS.NMDS.png",
format="png", height=15, width=25, dpi=300, units="cm")


# plot time only
P9 = (
ggplot(data=x_time)
+ geom_point(mapping=aes(x="X", y="Y", color="Timepoint", fill="Timepoint"), size=3, show_legend=None)
+ stat_ellipse(	mapping=aes(x="X", y="Y", color="Timepoint"), size=0.6, show_legend=None,
				type="t", level=0.95, segments=51)
+ ggtitle("")
+ xlab("NMDS1")
+ ylab("NMDS2")
+ scale_fill_manual(name="", values=TimeColors)
+ scale_color_manual(name="", values=TimeColors, guide=False)
+ scale_size_continuous(guide=False)
+ theme(plot_title=element_text(family="sans", face="bold",
					colour='#261A1A', size=12),
	axis_title_x=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_title_y=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_text_x=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	axis_text_y=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	legend_text=element_text(family="sans", face="plain", colour='#000000', size=10),
	legend_title=element_text(family="sans", face="bold", colour='#000000', size=10),
	panel_background=element_rect(fill="white", color="white"),
	panel_grid=element_line(color="#E6E6E6"),
	aspect_ratio=0.3,
	panel_spacing=0.5)
)

P9.save(filename="/binfl/lv70694/schmat/chicken/diversity_taxa/RESULTS.NMDS.timepoints.png",
format="png", height=15, width=25, dpi=300, units="cm")


# plot treat only
P10 = (
ggplot(data=x_treat)
+ geom_point(mapping=aes(x="X", y="Y", color="Treatment", fill="Treatment"), size=3, show_legend=None)
+ stat_ellipse(	mapping=aes(x="X", y="Y", color="Treatment"), size=0.6, show_legend=None,
				type="t", level=0.95, segments=51)
+ ggtitle("")
+ xlab("NMDS1")
+ ylab("NMDS2")
+ scale_fill_manual(name="", values=TreatColors)
+ scale_color_manual(name="", values=TreatColors, guide=False)
+ scale_size_continuous(guide=False)
+ theme(plot_title=element_text(family="sans", face="bold",
					colour='#261A1A', size=12),
	axis_title_x=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_title_y=element_text(	family="sans", face="bold",
					colour='#261A1A', size=9),
	axis_text_x=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	axis_text_y=element_text(family="sans", face="plain", colour='#261A1A', size=11),
	legend_text=element_text(family="sans", face="plain", colour='#000000', size=10),
	legend_title=element_text(family="sans", face="bold", colour='#000000', size=10),
	panel_background=element_rect(fill="white", color="white"),
	panel_grid=element_line(color="#E6E6E6"),
	aspect_ratio=0.3,
	panel_spacing=0.5)
)

P10.save(filename="/binfl/lv70694/schmat/chicken/diversity_taxa/RESULTS.NMDS.treatments.png",
format="png", height=15, width=25, dpi=300, units="cm")
