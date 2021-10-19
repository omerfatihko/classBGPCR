import pandas as pd
import time
start_time = time.time()
pd.set_option("display.max_rows", None, "display.max_columns", None)


#read the file
allvs23comp = pd.read_csv("/cta/users/ofkonar/work/results/csvs/RAMPall_vs_2-3.csv")

#count number of speciality occurrences for each position and add it to the result df
columnnames = list(allvs23comp.columns)
#get only blosum score columns
blosums = [x for x in columnnames if "blosum" in x]
blosumdf = allvs23comp[blosums]
#we turn each row to a list, count number of blosum scores and add it to the df
allvs23comp = pd.concat([allvs23comp, blosumdf.apply(lambda row: len([x for x in list(row) if x != "-"]), axis = 1)], axis = 1)
allvs23comp.rename(columns = {0 : "count"}, inplace = True)
print(allvs23comp.head())
allvs23comp.to_csv("/cta/users/ofkonar/work/results/csvs/RAMPall_vs_2-3_binary_comparisons.csv") #, index = False


#read the file
allvs23comp = pd.read_csv("/cta/users/ofkonar/work/results/csvs/RAMPall_vs_2-3_binary_comparisons.csv")
#select the residues you want to use
selected = allvs23comp.loc[allvs23comp["count"]>16, "calrl_seq"].tolist()
print(selected)
#read the score file 
score = pd.read_csv("/cta/users/ofkonar/work/structures/all_structures/6e3y.pdb.cscore.csv")
#get the subset of interactions that involves at least one receptor residue
receptorresidues = score.loc[(score["chain_1"] == "R") | (score["chain_2"] == "R")]
#filter the selected residues
filter1 = receptorresidues["residue_1"].isin(selected)
filter2 = receptorresidues["residue_2"].isin(selected)
selectedinteractions = receptorresidues[filter1 | filter2]

selectedinteractions.to_csv("/cta/users/ofkonar/work/results/csvs/RAMPall_vs_2-3_selected_interactions.csv", index = False)


print("My program took", time.time() - start_time, "seconds to run")

