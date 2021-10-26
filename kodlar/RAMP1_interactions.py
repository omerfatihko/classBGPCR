from functions import my_functions as mf
import pandas as pd
import time
from os import access, listdir
from os.path import isfile, join
pd.set_option("display.max_rows", None, "display.max_columns", None)

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
scoredf1 = pd.read_csv("/cta/users/ofkonar/work/structures/all_structures/6e3y.pdb.cscore.csv")
filter1 = scoredf1.loc[(scoredf1["chain_1"] == "RAMP1") & (scoredf1["chain_2"] == "R")]
reslist1 = filter1["residue_2"].tolist()
filter2 = scoredf1.loc[(scoredf1["chain_1"] == "R") & (scoredf1["chain_2"] == "RAMP1")]
reslist2 = filter2["residue_1"].tolist()
collectedlist1 = set(reslist1 +reslist2)
print(len(collectedlist1))
#get the subset of residues that interact with RAMP1
selected1  = globalconsensus.loc[globalconsensus["calrl_seq"].isin(collectedlist1)]
print(selected1.shape)
ramp1indexes = globalconsensus.index[globalconsensus["calrl_seq"].isin(collectedlist1)].tolist()
print(ramp1indexes)
#selected1.to_csv(csvpath + "selected1.csv", index=False)

#read 6UUN score file to get residues interacting with RAMP2
scoredf2 = pd.read_csv("/cta/users/ofkonar/work/structures/all_structures/6uun.pdb.cscore.csv")
filter3 = scoredf2.loc[(scoredf2["chain_1"] == "RAMP2") & (scoredf2["chain_2"] == "R")]
reslist3 = filter3["residue_2"].tolist()
filter4 = scoredf2.loc[(scoredf2["chain_1"] == "R") & (scoredf2["chain_2"] == "RAMP2")]
reslist4 = filter4["residue_1"].tolist()
collectedlist2 = set(reslist3 +reslist4)

#get the subset of residues that interact with RAMP1
selected2  = globalconsensus.loc[globalconsensus["calrl_seq"].isin(collectedlist2)]
#selected2.to_csv(csvpath + "selected2.csv", index=False)

#read 6UUS and 6UVA score files to get residues interacting with RAMP3
scoredf3 = pd.read_csv("/cta/users/ofkonar/work/structures/all_structures/6uus.pdb.cscore.csv")
scoredf4 = pd.read_csv("/cta/users/ofkonar/work/structures/all_structures/6uva.pdb.cscore.csv")
filter5 = scoredf3.loc[(scoredf3["chain_1"]=="RAMP3") & (scoredf3["chain_2"]=="R")]["residue_2"].tolist()
filter6 = scoredf3.loc[(scoredf3["chain_1"]=="R")&(scoredf3["chain_2"]=="RAMP3")]["residue_1"].tolist()
collected6uus = sorted(list(set(filter5+filter6)))
filter7 = scoredf4.loc[(scoredf4["chain_1"]=="RAMP3")&(scoredf4["chain_2"]=="R")]["residue_2"].tolist()
filter8 = scoredf4.loc[(scoredf4["chain_1"]=="R")&(scoredf4["chain_2"]=="RAMP3")]["residue_1"].tolist()
collected6uva = sorted(list(set(filter7+filter8)))
collectedlist3 = set(collected6uva+collected6uus)

#get the subset of residues that interact with RAMP1
selected3  = globalconsensus.loc[globalconsensus["calrl_seq"].isin(collectedlist3)]
selected3.to_csv(csvpath + "selected3.csv", index=False)

print(collectedlist3)
print(collectedlist2)
print(collectedlist1)
allofthem = (collectedlist1.intersection(collectedlist2)).intersection(collectedlist3)
print("intersection of all three: ", allofthem)
print("------------------------------------------------------------------------------------------")
only1 = (collectedlist1.difference(collectedlist2)).difference(collectedlist3)
print("only ramp1: ", only1)
print("------------------------------------------------------------------------------------------")
only1and2 = (collectedlist1.intersection(collectedlist2)).difference(collectedlist3)
print("only 1 and 2 not 3", only1and2)
print("------------------------------------------------------------------------------------------")
only1and3 = (collectedlist1.intersection(collectedlist3)).difference(collectedlist2)
print("only 1 and 3 not 2: ", only1and3)
print("------------------------------------------------------------------------------------------")
only2 = (collectedlist2.difference(collectedlist1)).difference(collectedlist3)
print("only 2:", only2)
print("------------------------------------------------------------------------------------------")
only2and3 = (collectedlist2.intersection(collectedlist3)).difference(collectedlist1)
print("only 2 and 3 not 1:", only2and3)
print("------------------------------------------------------------------------------------------")
only3 = (collectedlist3.difference(collectedlist2)).difference(collectedlist1)
print("only 3: ", only3)
print("------------------------------------------------------------------------------------------")

#check whether they share a common aminoacid for said residue positions
allclassB = ["calcr", "calrl", "crfr1", "crfr2", "gcgr", "ghrhr", "gipr", "glp1r", "glp2r", "pacr", "pth1r", "pth2r", "sctr", "vipr1","vipr2"]
ramp1interacting = ["calcr", "calrl", "gcgr", "gipr", "glp1r", "glp2r", "pacr", "pth1r", "pth2r", "sctr", "vipr2"]
ramp1notinteracting = ["crfr1", "ghrhr", "vipr1"]
onlyramp3interacting = ["crfr2"]
fastapath = "/cta/users/ofkonar/work/resources/class_B1/canonical/"
onlyfiles = [f for f in listdir(fastapath) if isfile(join(fastapath, f))]
orthologfiles = sorted([f for f in onlyfiles if mf.FilterSeq(f, allclassB) if "orthologs" in f])
fastadict = {}
for i in orthologfiles:
    prtname = i.split("_")[0]
    fastadict[prtname] = mf.read_fasta(fastapath+i)

ramp1interactingfastas = [fastadict.get(key) for key in ramp1interacting]
flattenedramp1fastalist = [item for sublist in ramp1interactingfastas for item in sublist]
ramp1interactingdf = mf.fasta_to_dataframe(flattenedramp1fastalist)
ramp1consensus = mf.consensus(ramp1interactingdf, 0.9)
print(ramp1consensus.iloc[ramp1indexes])
print("------------------------------------------------------------------------------------------")
a = mf.residuecontent(ramp1interactingdf)
#print(a)
print(a.iloc[ramp1indexes])


print("My program took", time.time() - start_time, "seconds to run")