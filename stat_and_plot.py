#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
import glob

datapath = sys.argv[1]

meta = pd.read_csv(str(datapath) + "/param_inputs_clean.txt", sep = "\s", names = ["map", "SN", "initial_ratio"])
meta = meta.reset_index(drop=False)
meta.columns = ["run", "map", "SN", "initial_ratio"]
meta["run"] = meta["run"] + 1

meta["map"] = meta['map'].str.split('=', 1, expand=True)[1]
meta["SN"] = meta['SN'].str.split('=', 1, expand=True)[1].astype(float)
meta["initial_ratio"] = meta['initial_ratio'].str.split('=', 1, expand=True)[1].astype(float)

empty = pd.DataFrame(columns = ['Year', 'PopulationSize', 'Total_palmata', 'Total_cervicornis',
       'Total_prolifera', 'Ratio_palmata_to_cervicornis',
       'cervicornis_ancestry_in_palmata', 'palmata_ancestry_in_cervicornis',
       'palmata_ancestry_in_prolifera', 'cervicornis_ancestry_in_prolifera', "run"])

for i in range(1,len(glob.glob(str(datapath) + "/sim_pop_stats_per_year_run*"))+1):
    tmp = pd.read_csv(str(datapath) + "/sim_pop_stats_per_year_run" + str(i) + ".csv")
    tmp["run"] = i
    empty = pd.concat([empty, tmp]).reset_index(drop=True)

empty = empty.merge(meta, how = "left", on = "run")

unique_starting = empty.groupby(['SN','initial_ratio']).size().reset_index()
unique_starting = unique_starting[["SN", "initial_ratio"]]

def f_test(x, y):
    x = np.array(x)
    y = np.array(y)
    f = np.var(x, ddof=1)/np.var(y, ddof=1) #calculate F test statistic 
    dfn = x.size-1 #define degrees of freedom numerator 
    dfd = y.size-1 #define degrees of freedom denominator 
    p = 1-scipy.stats.f.cdf(f, dfn, dfd) #find p-value of F test statistic 
    return f, p

simple = empty[empty["map"] == "Simple"].reset_index(drop=True)
caye = empty[empty["map"] == "Caye"].reset_index(drop=True)

def pairwise_f_and_ttest(frame, outprefix):
    unique_comparisons = pd.DataFrame(columns = ["x_SN_IR", "y_SN_IR"])
    for i in range(len(unique_starting)):
        x_SN = unique_starting["SN"][i]
        x_IR = unique_starting["initial_ratio"][i]
        comparisons = list(range(len(unique_starting)))
        comparisons.remove(i)
        for compare in comparisons:
            y_SN = unique_starting["SN"][compare]
            y_IR = unique_starting["initial_ratio"][compare]
            unique_comparisons = pd.concat([unique_comparisons,pd.DataFrame([[str(x_SN)+"_"+str(x_IR), str(y_SN)+"_"+str(y_IR)]], columns = ["x_SN_IR", "y_SN_IR"])]).reset_index(drop=True)
    unique_comparisons.values.sort()
    unique_comparisons = unique_comparisons.drop_duplicates(["x_SN_IR", "y_SN_IR"]).reset_index(drop=True)
    stat_tests = pd.DataFrame(columns = ["generation", "SN_1", "IR_1", "SN_2", "IR_2", "fstat", "fstat_pval", "tstat", "tstat_pval"])
    for gen in [200, 10000]:
        for i in range(len(unique_comparisons)):
            x_SN = float(unique_comparisons["x_SN_IR"][i].split('_')[0])
            x_IR = float(unique_comparisons["x_SN_IR"][i].split('_')[1])
            y_SN = float(unique_comparisons["y_SN_IR"][i].split('_')[0])
            y_IR = float(unique_comparisons["y_SN_IR"][i].split('_')[1])
            x_chains = np.array(frame["palmata_ancestry_in_cervicornis"][frame["initial_ratio"] == x_IR][frame["SN"] == x_SN][frame["Year"] == gen])
            y_chains = np.array(frame["palmata_ancestry_in_cervicornis"][frame["initial_ratio"] == y_IR][frame["SN"] == y_SN][frame["Year"] == gen])
            t_result = scipy.stats.ttest_ind(x_chains, y_chains, equal_var=False)
            f_result = f_test(x_chains, y_chains)
            stat_tests = pd.concat([stat_tests, pd.DataFrame([[gen, x_SN, x_IR, y_SN, y_IR, f_result[0], f_result[1], t_result[0], t_result[1]]], columns = ["generation", "SN_1", "IR_1", "SN_2", "IR_2", "fstat", "fstat_pval", "tstat", "tstat_pval"])]).reset_index(drop=True)
    stat_tests.to_csv(str(datapath) + "/" + str(outprefix) + "_map_stats.csv")

pairwise_f_and_ttest(simple, "simple")
pairwise_f_and_ttest(caye, "caye")

def plot_runs(df, outprefix):
    for gen in [200, 10000]:
        for i in range(len(unique_starting)):
            sns.set_theme(style="ticks")
            SN = unique_starting["SN"][i]
            IR = unique_starting["initial_ratio"][i]
            plt.figure()
            sns.set_theme(style="ticks")
            # Define the palette as a list to specify exact values
            palette = sns.color_palette("rocket_r")
            max_y = max(df["cervicornis_ancestry_in_palmata"][df["Year"] <= gen].tolist() + df["palmata_ancestry_in_cervicornis"][df["Year"] <= gen].tolist())
            tmp = df[df["initial_ratio"] == IR][df["SN"] == SN][df["Year"] <= gen].reset_index(drop=True)
            tmp["Ancestry"] = tmp["palmata_ancestry_in_cervicornis"]
            # Plot the lines on two facets
            sns.lineplot(
                data=tmp,
                x="Year", y="Ancestry",
                hue="run", palette=palette
            )
            tmp1 = df[df["initial_ratio"] == IR][df["SN"] == SN][df["Year"] <= gen].reset_index(drop=True)
            tmp1["Ancestry"] = tmp1["cervicornis_ancestry_in_palmata"]
            palette = sns.color_palette("mako_r")
            sns.lineplot(
                data=tmp1,
                x="Year", y="Ancestry",
                hue="run", palette=palette
            )

            plt.ylim(0, max_y + 0.1*max_y)
            plt.legend([],[], frameon=False)
            plt.title('SN = ' + str(SN) + '; Initial Ratio = ' + str(IR), y=1.0, pad=-14)
            plt.tight_layout()
            plt.savefig(str(datapath) + "/" + str(outprefix) + "_YR" + str(gen) + "_SN" + str(SN) + "_IR" + str(IR) + ".png", dpi = 300)

plot_runs(simple, "simple")
plot_runs(caye, "caye")
