from numpy import empty
from functions import my_functions as mf
import pandas as pd
import time
from os import access, listdir
from os.path import isfile, join
pd.set_option("display.max_rows", None, "display.max_columns", None)

start_time = time.time()

#read the files that would be used most of the analyses, define the directories that will be used
B1consensus = pd.read_csv("/cta/users/ofkonar/work/results/csvs/class_B1_all_consensus_seqs.csv")
fastadir = "/cta/users/ofkonar/work/resources/class_B1/canonical/"
onlyfiles = [f for f in listdir(fastadir) if isfile(join(fastadir, f))] #files in fasta directory

#define the subfamilies
gcgrsubfamily = {"gcgr", "gipr", "glp1r", "glp2r"}
sctrsubfamily = {"sctr", "pacr", "ghrhr", "vipr1", "vipr2"}
calcrsubfamily = {"calcr", "calrl", "crfr1", "crfr2"}
pthfamily = {"pth1r", "pth2r"}

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

print("My program took", time.time() - start_time, "seconds to run")