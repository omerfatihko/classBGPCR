from functions import my_functions as mf
import pandas as pd
import time
from os import access, listdir
from os.path import isfile, join
start_time = time.time()

#first we need to gather all local consensus with global consensus so we know which residue from one GPCR corresponds
#to which residue from another GPCR.

#read the global consensus 
globalconsensus = pd.read_csv("/cta/users/ofkonar/work/results/csvs/class_B1_global_consensus.csv", index_col = 0)
#Binary comparison files already includes the GPCR sequences aligned with the global consensus, so we will loop over them and take
#appropriate columns (3,4,5, and 6th columns)
GPCRlist = ["calcr", "calrl", "crfr1", "crfr2", "gcgr", "ghrhr", "gipr", "glp1r", "glp2r", "pacr", "pth1r", "pth2r", "sctr", "vipr1", "vipr2"]
csvpath = "/cta/users/ofkonar/work/results/csvs/"
onlyfiles = [f for f in listdir(csvpath) if isfile(join(csvpath, f))]
bincompfiles = sorted([i for i in onlyfiles if mf.FilterSeq(i, GPCRlist) if "binary" in i])
for f in bincompfiles:
    tempdf = pd.read_csv(csvpath + f)
    localconsensus = tempdf.iloc[:,[3,4,5,6]]
    globalconsensus = pd.concat([globalconsensus, localconsensus], axis=1)
#globalconsensus.to_csv(csvpath + "class_B1_all_consensus_seqs.csv", index=False)

#read 6E3Y score file to get residues interacting with RAMP1
scoredf = pd.read_csv("/cta/users/ofkonar/work/structures/all_structures/6e3y.pdb.cscore.csv")
filter1 = scoredf.loc[(scoredf["chain_1"] == "RAMP1") & (scoredf["chain_2"] == "R")]
reslist1 = filter1["residue_2"].tolist()
filter2 = scoredf.loc[(scoredf["chain_1"] == "R") & (scoredf["chain_2"] == "RAMP1")]
reslist2 = filter2["residue_1"].tolist()
reslist = sorted(list(set(reslist1 +reslist2)))

#get the subset of residues that interact with RAMP1
selected  = globalconsensus.loc[globalconsensus["calrl_seq"].isin(reslist)]
#selected.to_csv(csvpath + "selected.csv", index=False)

ramp1interacting = ["calcr", "calrl", "gcgr", "gipr", "glp1r", "glp2r", "pacr", "pth1r", "pth2r", "sctr", "vipr2"]
ramp1notinteracting = ["crfr1", "ghrhr", "vipr1"]
onlyramp3interacting = ["crfr2"]

ramp1orthologs = sorted([f for f in onlyfiles if mf.FilterSeq(f, ramp1interacting) if "orthologs" in f])
print(ramp1orthologs)
print("My program took", time.time() - start_time, "seconds to run")