from numpy import empty
from functions import my_functions as mf
import pandas as pd
import time
from os import access, listdir, path
from os.path import isfile, join
pd.set_option("display.max_rows", None, "display.max_columns", None)

start_time = time.time()

#read the files that would be used most of the analyses, define the directories that will be used
B1consensus = pd.read_csv("/cta/users/ofkonar/work/results/csvs/class_B1_all_consensus_seqs.csv")
fastadir = "/cta/users/ofkonar/work/resources/class_B1/canonical/"
#csvdir = "/cta/users/ofkonar/work/results/csvs/"
onlyfiles = [f for f in listdir(fastadir) if isfile(join(fastadir, f))] #files in fasta directory

#define the subfamilies
calcrsubfamily = {"calcr", "calrl", "crfr1", "crfr2"}
pthsubfamily = {"pth1r", "pth2r"}
sctrsubfamily = {"sctr", "pacr", "ghrhr", "vipr1", "vipr2"}
gcgrsubfamily = {"gcgr", "gipr", "glp1r", "glp2r"}

#################################################################################################
"""clades are collected as fasta files
#to define, we have to construct tree of class B
#clades are already collected from previously constructed trees
#collect the clades into one list
cladefiles = [f for f in onlyfiles if "_clade.txt" in f]
#collect the fasta files for each clade, write them down
vipr1file = [k for k in cladefiles if "vipr1" in k]
for i in vipr1file:
    prtname = i.split("_")[0]
    targetfastadir = [j for j in onlyfiles if prtname in j if "_blast.fasta" in j]
    #read the txt file and convert it to a list
    with open(fastadir + i) as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]
    #filter the clade from the target fasta
    cladefasta = mf.filter_fasta(lines, fastadir + targetfastadir[0])
    #write the fasta list down
    mf.write_fasta(fastadir + prtname + "_clade.fasta", cladefasta)
"""
#################################################################################################
""" clades are combined and aligned
#combine clade fastas into one file and align them
cladefastas = [f for f in onlyfiles if "_07" in f if "clstr" not in f if "B1" not in f]
allcladefastaslist = []
for i in cladefastas:
    temp = mf.read_fasta(fastadir + i)
    allcladefastaslist = allcladefastaslist + temp
print(len(allcladefastaslist))
setfasta = set(allcladefastaslist)
print(len(setfasta))
mf.write_fasta(fastadir + "Class_B1_cdhit_07.fasta", allcladefastaslist)
"""
#################################################################################################
"""
#there were 2 duplicate sequences in the alignment, eliminate them, since we are interested in the
#general shape of the tree, that won't make much of a change on it.
#read the alignment file
classB1fastas = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/Class_B1_cdhit_07.fasta")
classB1set = list(set(classB1fastas))
print(len(classB1set))
print(len(classB1fastas))
mf.write_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/Class_B1_cdhit_07_withoutgroup.fasta", classB1set)
"""
#################################################################################################
#calcr subfamily versus pth subfamily
#this part is quarantined because below technique is better (compare one by one rather than comparing consensus of whole family)
"""
#gather calcr subfamily sequences and form a subfamily consensus
calcrsubfastalist = []
for i in calcrsubfamily:
    targetfasta = [f for f in onlyfiles if i in f if "msa_pruned" in f]
    targetfastalist = mf.read_fasta(fastadir + targetfasta[0])
    calcrsubfastalist = calcrsubfastalist + targetfastalist
#get the consensus of the subfamily
calcrsubdf = mf.fasta_to_dataframe(calcrsubfastalist)
calcrsubcon = mf.consensus(calcrsubdf, 0.9)
#gather pth subfamily sequences and form a subfamily consensus
pthsubfastalist = []
for i in pthsubfamily:
    targetfasta = [f for f in onlyfiles if i in f if "msa_pruned" in f]
    targetfastalist = mf.read_fasta(fastadir + targetfasta[0])
    pthsubfastalist = pthsubfastalist + targetfastalist
#get the consensus of the subfamily
pthsubdf = mf.fasta_to_dataframe(pthsubfastalist)
pthsubcon = mf.consensus(pthsubdf, 0.9)

#compare calcr subfamily to pth subfamily
calcrsubspecific = mf.special_residues(calcrsubcon, pthsubcon)
#print(calcrsubspecific)

#find the residues on the receptors
calcrsubspecidx = calcrsubspecific.loc[calcrsubspecific["residue"] != "-"].index.tolist()
targetresidues = B1consensus.iloc[calcrsubspecidx]

#whole consensus vs whole consensus
targetresidues.to_csv("/cta/users/ofkonar/work/tae/birbakalim.csv", index = False)
"""
"""
#BURAYI FOR LOOP İLE TEKRAR YAPMAK LAZIM, HER SEFERİNDE TEK TEK RESEPTÖR İSMİ Mİ YAZACAĞIM!
#compare calcr subfamily members between themselves and keep the residues that agree for all
#compare calcr subfamily members and pth subfamily members one by one and keep residues that disagree for all
#get relevant part of the B1 consensus dataframe
colnames = B1consensus.columns.tolist()
relevantnames = [i for i in colnames if mf.FilterSeq(i, list(calcrsubfamily) + list(pthsubfamily)) == i]
relevantdf = B1consensus[relevantnames]
#first get calcr subfamily residues where all are conserved for each member
calcrC = relevantdf.loc[(relevantdf["calcr_status"]=="C")&(relevantdf["calrl_status"]=="C")&(relevantdf["crfr1_status"]=="C")&(relevantdf["crfr2_status"]=="C")]
#next, get residues where calcr subfamily agrees
calcragree = calcrC.loc[(calcrC["calcr_residue"]==calcrC["calrl_residue"])&(calcrC["calcr_residue"]==calcrC["crfr1_residue"])&(calcrC["calcr_residue"]==calcrC["crfr2_residue"])]
#then get residues where calcr subfamily disagrees with pth subfamily
calcrdifpth = calcragree.loc[(calcragree["calcr_residue"]!=calcragree["pth1r_residue"])&(calcragree["calcr_residue"]!=calcragree["pth2r_residue"])]
#print(calcrdifpth)

#calcr family vs others
#conserved for all calcr family
#get calcr rows that are consensus (rows we care)
calcrstatuscols = [i for i in colnames if "status" in i if mf.FilterSeq(i, list(calcrsubfamily)) == i]
calcrrescols = [i for i in colnames if "residue" in i if mf.FilterSeq(i, list(calcrsubfamily)) == i]
restrescols = [i for i in colnames if "_residue" in i if mf.FilterSeq(i, list(calcrsubfamily)) != i]
#get rows calcr subfamily residues are all conserved "C"
subdf = B1consensus.copy()
for i in calcrstatuscols:
    subdf = subdf.loc[subdf[i] == "C"]
#get rows calcr subfamily residues agree with each other
#sliding windows of pairs compared with each other
for i in range(len(calcrrescols)-1):
    pair = calcrrescols[i:i+2]
    subdf = subdf.loc[subdf[pair[0]] == subdf[pair[1]]]
#get rows calcr subfamily is different from the rest
for i in restrescols:
    subdf = subdf.loc[subdf[calcrrescols[0]] !=subdf[i]]
print(subdf)
print("######################################################################################################################")
#calcrs vs crfrs
#conserved for both (and different than the rest) but different between the two
#conserved within calcr and calrl
group1 = B1consensus.loc[(B1consensus["calcr_status"] == "C") & (B1consensus["calrl_status"] == "C") & (B1consensus["calcr_residue"] == B1consensus["calrl_residue"])]
group1 = group1.loc[(group1["crfr1_status"] == "C") & (group1["crfr2_status"] == "C") & (group1["crfr1_residue"] == group1["crfr2_residue"])]
for i in restrescols:
    group1 = group1.loc[group1["calcr_residue"] != group1[i]]
    group1 = group1.loc[group1["crfr1_residue"] != group1[i]]
group1 = group1.loc[group1["calcr_residue"] != group1["crfr1_residue"]]
print(group1)
print("######################################################################################################################")
#conserved for calcrs but not for crfrs
group2 = B1consensus.loc[(B1consensus["calcr_status"] == "C") & (B1consensus["calrl_status"] == "C") & (B1consensus["calcr_residue"] == B1consensus["calrl_residue"]) & (B1consensus["crfr1_status"] == "NC") & (B1consensus["crfr2_status"] == "NC")]
print(group2)
print("######################################################################################################################")
"""
#define the class B1 family tree
#        _______n0______
#       /               \
#    n1            ______n2________
#    /\           /                \
#glp2r n3       n4              ____n5___
#      /\      / \             /         \
#  glp1r n6 sctr  n7         n8        ___n9___
#        /\       /\         /\       /        \
#    gipr gcgr pacr n10 pth1r pth2r n11         n12
#                   /\              /\          /\
#               vipr2 n13       crfr2 crfr1 calcr calrl
#                     /\  
#                 ghrhr vipr1

ghrhr = mf.Node("ghrhr")
vipr1 = mf.Node("vipr1")
vipr2 = mf.Node("vipr2")
n13 = mf.Node("n13", left = ghrhr, right = vipr1)
crfr2 = mf.Node("crfr2")
crfr1 = mf.Node("crfr1")
calcr = mf.Node("calcr")
calrl = mf.Node("calrl")
gipr = mf.Node("gipr")
gcgr = mf.Node("gcgr")
pacr = mf.Node("pacr")
n10 = mf.Node("n10", left = vipr2, right = n13)
pth1r = mf.Node("pth1r")
pth2r = mf.Node("pth2r")
n11 = mf.Node("n11", left = crfr2, right = crfr1)
n12 = mf.Node("n12", left = calcr, right = calrl)
glp1r = mf.Node("glp1r")
n6 = mf.Node("n6", left = gipr, right = gcgr)
sctr = mf.Node("sctr")
n7 = mf.Node("n7", left = pacr, right = n10)
n8 = mf.Node("n8", left = pth1r, right = pth2r)
n9 = mf.Node("n9", left = n11, right = n12)
glp2r = mf.Node("glp2r")
n3 = mf.Node("n3", left = glp1r, right = n6)
n4 = mf.Node("n4", left = sctr, right = n7)
n5 = mf.Node("n5", left = n8, right = n9)
n1 = mf.Node("n1", left = glp2r, right = n3)
n2 = mf.Node("n2", left = n4, right = n5)
n0 = mf.Node("n0", left = n1, right = n2)

leaflist = []
mf.GetLeafList(n0, leaflist)
#print(leaflist)
#print("######################################################################################################################")

colnames = B1consensus.columns.tolist()
for i in leaflist:
    targetcolnames = [j for j in colnames if i.name in j]
    data = B1consensus.loc[:,targetcolnames]
    i.data = data
    #print(i.data)
#print("######################################################################################################################")

leafdictpopulated = {}
mf.GetLeafData(n9, leafdictpopulated)
print(leafdictpopulated)
print("######################################################################################################################")

print(mf.innercheck(leafdictpopulated))
a = mf.innercheck(leafdictpopulated)
print(B1consensus.iloc[a])

print("My program took", time.time() - start_time, "seconds to run")