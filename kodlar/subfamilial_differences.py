#import re
#from tracemalloc import DomainFilter
#from numpy import empty, inner
from functions import my_functions as mf
import pandas as pd
import time
from os import listdir #X_OK, access, ctermid, path
from os.path import isfile, join
#import statistics as stat
#from statistics import mean
#from kodlar.functions.my_functions import inner_check
#pd.set_option("display.max_rows", None, "display.max_columns", None)

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

#returns the list of leaf node objects
leaflist = mf.get_leaf_list(n0)
#fill the leaves with the data (dataframe)
colnames = B1consensus.columns.tolist()
for i in leaflist:
    targetcolnames = [j for j in colnames if i.name in j]
    data = B1consensus.loc[:,targetcolnames]
    i.data = data

#returns list of all node objects
allnodes = mf.get_node_list(n0)

for i in allnodes:
    #fill all nodes with inner check results
    nodedata = mf.get_leaf_data(i) # get the data of all leaves of the node
    nodeinnercheck = mf.inner_check(nodedata) #get inner check indexes
    nodemetadata = {"inner_check" : nodeinnercheck} #form a metadata dict
    i.metadata = nodemetadata #assign metadata dict to the node

for i in allnodes:    
    #fill all nodes with emergence, divergence, and difference check results (various comparisons of sister nodes)
    left = i.left #left sister
    right = i.right #right sister
    leftdata = mf.get_leaf_data(left) #get leaf data of the left sister
    rightdata = mf.get_leaf_data(right) #get leaf data of the right sister
    selfleaves = mf.get_leaf_list(i) #get leaf list of the major node
    selfleavesdata = mf.get_leaf_data(i) #get leaf data of the major node
    rest = [j for j in leaflist if j not in selfleaves] #get leaves that do not belong to the target node
    restdata = {}
    for j in rest:
        restdata.update(mf.get_leaf_data(j)) #get leaf data of rest of the leaves (do not belong to the major node)
    
    if left and right: #node should have children (not leaf), no need to check for metadata anymore, there is some for sure
        lefttorightdiv = mf.divergence(leftdata, rightdata) #divergence of the left sister
        leftorightem = mf.emergence(leftdata, rightdata) #emergence of the left sister
        lefttorightdif = mf.outer_check(leftdata, rightdata) #difference of the left sister
        left.metadata["divergence"] = lefttorightdiv
        left.metadata["emergence"] = leftorightem
        left.metadata["difference"] = lefttorightdif
        righttoleftdiv = mf.divergence(rightdata, leftdata) #divergence of the right sister
        righttoleftem = mf.emergence(rightdata, leftdata) #emergence of the right sister
        righttoleftdif = mf.outer_check(rightdata, leftdata) #difference of the right sister 
        right.metadata["divergence"] = righttoleftdiv
        right.metadata["emergence"] = righttoleftem
        right.metadata["difference"] = righttoleftdif
        if restdata: #if restdata is not empty (if i is not root)
            leftoutercheck = mf.outer_check(leftdata, restdata) #outer check of the left sister
            left.metadata["outer_check"] = leftoutercheck
            rightoutercheck = mf.outer_check(rightdata, restdata) #outer check of the right sister
            right.metadata["outer_check"] = rightoutercheck

#add Ogun Hoca's outer check idea (closest that is not sister as outgroup)
for i in allnodes:
    left = i.left
    right = i.right
    if left and right: #first, a node should have 2 children
        #for left
        leftsleft = left.left
        leftsright = left.right
        if leftsleft or leftsright: # if left node has any children
            #instead of rest of the tree becoming the control for outer check, aunt node becomes the outgroup (sister of the ancestor node)
            outgroupdata = mf.get_leaf_data(right)
            if leftsleft:
                targetdata = mf.get_leaf_data(leftsleft)
                smalloutercheck = mf.outer_check(targetdata, outgroupdata)
                leftsleft.metadata["Ogun's_outer_check"] = smalloutercheck
            if leftsright:
                targetdata = mf.get_leaf_data(leftsright)
                smalloutercheck = mf.outer_check(targetdata, outgroupdata)
                leftsright.metadata["Ogun's_outer_check"] = smalloutercheck
        #for right
        rightsleft = right.left
        rightsright = right.right
        if rightsleft or rightsright:
            outgroupdata = mf.get_leaf_data(left)
            if rightsleft:
                targetdata = mf.get_leaf_data(rightsleft)
                smalloutercheck = mf.outer_check(targetdata, outgroupdata)
                rightsleft.metadata["Ogun's_outer_check"] = smalloutercheck
            if rightsright:
                targetdata =mf.get_leaf_data(rightsright)
                smalloutercheck = mf.outer_check(targetdata, outgroupdata)
                rightsright.metadata["Ogun's_outer_check"] = smalloutercheck

#for i in allnodes:
#    print(i.name)
#    metadata = i.metadata
#    for key in metadata:
#        print(key)
#        print(metadata[key])
#    print("\n")

for i in allnodes:
    print(i.name)
    #emergence
    #without outer check
    if "emergence" in i.metadata:
        nakedemergence = [j for j in i.metadata["emergence"] if j in i.metadata["inner_check"]]
        print("emergence without outer check: ")
        print(nakedemergence)
    #with outer check
    if "outer_check" in i.metadata:
        temp1 = [i.metadata["inner_check"], i.metadata["emergence"], i.metadata["outer_check"]]
        dressedemergence = set(temp1[0]).intersection(*temp1)
        print("emergence with outer check: ")
        print(dressedemergence)
    #with Ogun's outer check
    if "Ogun's_outer_check" in i.metadata:
        temp11 = [i.metadata["inner_check"], i.metadata["emergence"], i.metadata["Ogun's_outer_check"]]
        ogunsemergence = set(temp11[0]).intersection(*temp11)
        print("emergence with Ogun's outer check: ")
        print(ogunsemergence)
        print("\n")
    #divergence
    #without outer check
    if "divergence" in i.metadata:
        nakeddivergence = [j for j in i.metadata["divergence"] if j in i.metadata["inner_check"]]
        print("divergence without outer check: ")
        print(nakeddivergence)
    #with outer check
    if "outer_check" in i.metadata:
        temp2 = [i.metadata["inner_check"], i.metadata["divergence"], i.metadata["outer_check"]]
        dresseddivergence = set(temp2[0]).intersection(*temp2)
        print("divergence with outer check: ")
        print(dresseddivergence)
    #with Ogun's outer check
    if "Ogun's_outer_check" in i.metadata:
        temp22 = [i.metadata["inner_check"], i.metadata["divergence"], i.metadata["Ogun's_outer_check"]]
        ogunsdivergence = set(temp22[0]).intersection(*temp22)
        print("divergence with Ogun's outer check: ")
        print(ogunsdivergence)
        print("\n")
    #difference
    #without outer check
    if "difference" in i.metadata:
        nakeddifference = [j for j in i.metadata["difference"] if j in i.metadata["inner_check"]]
        print("difference without outer check: ")
        print(nakeddifference)
    #with outer check
    if "outer_check" in i.metadata:
        temp3 = [i.metadata["inner_check"], i.metadata["difference"], i.metadata["outer_check"]]
        dresseddifference = set(temp3[0]).intersection(*temp3)
        print("difference with outer check: ")
        print(dresseddifference)
    #with Ogun's outer check
    if "Ogun's_outer_check" in i.metadata:
        temp33 = [i.metadata["inner_check"], i.metadata["difference"], i.metadata["Ogun's_outer_check"]]
        ogunsdifference = set(temp33[0]).intersection(*temp33)
        print("difference with Ogun's outer check: ")
        print(ogunsdifference)
        print("\n")
    print("######################################################################################################################")
    print("\n")
#statistic test, which domains (ECD, TMD, ICD) the residues come from?
#get domain info file
domaininfo = pd.read_csv("/cta/users/ofkonar/work/resources/class_B1_domains.csv", index_col = 0)
#ancestral trail of each receptor
glp2rlist = [glp2r, n1, n0]
glp1rlist = [glp1r, n3, n1, n0]
giprlist = [gipr, n6, n3, n1, n0]
gcgrlist = [gcgr, n6, n3, n1, n0]
sctrlist = [sctr, n4, n2, n0]
pacrlist = [pacr, n7, n4, n2, n0]
vipr2list = [vipr2, n10, n7, n4, n2, n0]
ghrhrlist = [ghrhr, n13, n10, n7, n4, n2, n0]
vipr1list = [vipr1, n13, n10, n7, n4, n2, n0]
pth1rlist = [pth1r, n8, n5, n2, n0]
pth2rlist = [pth2r, n8, n5, n2, n0]
crfr2list = [crfr2, n11, n9, n5, n2, n0]
crfr1list = [crfr1, n11, n9, n5, n2, n0]
calcrlist = [calcr, n12, n9, n5, n2, n0]
calrllist = [calrl, n12, n9, n5, n2, n0]
ancestrydict = {"glp2r" : glp2rlist, "glp1r" : glp1rlist, "gipr" : giprlist, "gcgr" : gcgrlist, "sctr" : sctrlist, "pacr" : pacrlist, "vipr2" : vipr2list, "ghrhr" : ghrhrlist, 
"vipr1" : vipr1list, "pth1r" : pth1rlist, "pth2r" : pth2rlist, "crfr2" : crfr2list, "crfr1" : crfr1list, "calcr" : calcrlist, "calrl" : calrllist}

#we will perform a chi square analysis to see whether a domain is favored at any level of the tree, for more on the chi square analysis
#visit https://www.statisticshowto.com/probability-and-statistics/chi-square/
for key in ancestrydict: #key is also the receptor's name
    ancestraltraillist = ancestrydict[key] #a leaf's ancestral trail (from leaf to the root)
    print(key)
    for node in ancestraltraillist:
        nodename = node.name
        print(nodename)
        if "Ogun's_outer_check" in node.metadata:
        #Difference with Ogun's outer check
            temp = [node.metadata["inner_check"], node.metadata["difference"], node.metadata["Ogun's_outer_check"]]
            ogunsdifference = list(set(temp[0]).intersection(*temp))
            ogunsdifference.sort()
            #get receptor length from domain info file
            reclength = domaininfo.loc["C-ter", key.upper()]
            #keep domain lengths in one list [N-ter, 7TM, C-ter]
            domainlenghts = []
            Nter = domaininfo.loc["N-ter", key.upper()]
            domainlenghts.append(Nter)
            TM7 = domaininfo.loc["TM7", key.upper()]
            corelength = TM7 - Nter
            domainlenghts.append(corelength)
            Cter = reclength - TM7
            domainlenghts.append(Cter)
            domainindexes = [] #where the domain ends
            domainindexes.append(domainlenghts[0])
            domainindexes.append(domainlenghts[0] + domainlenghts[1])
            domainindexes.append(domainlenghts[0] + domainlenghts[1] + domainlenghts[2])
            #turn the domain lengths into expected values(if we chose n residues how many would have
            # been from a given domain where n is number of residues we get from ogunsdifference)
            expected = [(x/reclength)*len(ogunsdifference) for x in domainlenghts]
            #print(expected)
           
            #now lets see what we observed
            observed = [0, 0, 0] #shares indexes with expected [N-ter, 7TM, C-ter]
            for i in ogunsdifference:
                residue = B1consensus.loc[i, key + "_seq"]
                if residue != "-":
                    resindex = int(residue[1:])
                    if resindex <= domainindexes[0]: #if residue index is within N-ter
                        observed[0] = observed[0] + 1
                    elif (domainindexes[0] < resindex) & (resindex <= domainindexes[1]):
                        observed[1] = observed[1] + 1
                    elif (domainindexes[1] < resindex) & (resindex <= domainindexes[2]):
                        observed[2] = observed[2] + 1
                    else:
                        print("receptor ", key, " has indexing issues")
                else:
                    observed[1] = observed[1] + 1
            #print(observed)
            if sum(observed) != 0:
                chiresults = mf.chisquare(observed, expected)
                chitable = pd.DataFrame([observed, expected],columns= ["N-ter", "7TM", "C-ter"])
                chitable.rename(index={0: "observed", 1: "expected"}, inplace=True)
                print(chitable.round(2))
                if chiresults[0] > 5.991:
                    print("p-value < 0.05 ", chiresults[0])
                else:
                    print("p-value is not significant")
            else:
                print("no observation")
        else:
            print("no Ogun's outer check")
                

    print("\n")
 print("######################################################################################################################")

#like in "/cta/users/ofkonar/work/tae/deneme2.py" use residue indexes to get a pymol script
#example script to select list of residues from a chain
"select deneme, chain A and resi 22+24+48+206"
"color white, chain A"
"hide (chain A)"
base_string = "select 11, chain R and resi "
for key in ancestrydict: #key is also the receptor's name
    ancestraltraillist = ancestrydict[key] #a leaf's ancestral trail (from leaf to the root)
    print(key)
    for node in ancestraltraillist:
        nodename = node.name
        print(nodename)
        if "Ogun's_outer_check" in node.metadata:
        #Difference with Ogun's outer check
            temp = [node.metadata["inner_check"], node.metadata["difference"], node.metadata["Ogun's_outer_check"]]
            ogunsdifference = list(set(temp[0]).intersection(*temp))
            ogunsdifference.sort()
            print(ogunsdifference)
            print(len(ogunsdifference), " residues")
            if ogunsdifference: #check if empty
                #get the selected residues for given receptor at each level (each ancestral node)
                selectedresidues = B1consensus.loc[ogunsdifference, key + "_seq"].tolist()
                selectedindexes = [x[1:] for x in selectedresidues if x != "-"]
                indexstring = ""
                for idx in selectedindexes:
                    indexstring += idx
                    indexstring += "+"
                indexstring = "select " + nodename + ", chain R and resi " + indexstring
                indexstring = indexstring[:-1]
                print(indexstring)
    print("\n")


#print(vars(n0))
#print(vars(gcgr))

#for i in allnodes:
#    print("node", i.name)
#    for key, value in i.metadata.items():
#        print(key)
#        print(value)
#    print("\n")

print("######################################################################################################################")
#subdf = B1consensus.loc[B1consensus["pacr_status"] == "C"]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["vipr2_status"] == "C"]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["ghrhr_status"] == "C"]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["vipr1_status"] == "C"]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["pacr_residue"] != subdf["sctr_residue"]]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["vipr2_residue"] != subdf["sctr_residue"]]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["ghrhr_residue"] != subdf["sctr_residue"]]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["vipr1_residue"] != subdf["sctr_residue"]]
#print(list(subdf.index))
#print(subdf)
#subdf = subdf.loc[subdf["crfr1_status"] == "C"]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["calcr_status"] == "C"]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["calrl_status"] == "C"]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["pacr_status"] == "NC"]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["vipr2_status"] == "NC"]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["ghrhr_status"] == "NC"]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["vipr1_status"] == "NC"]
#print(list(subdf.index))
#subdf = subdf.loc[subdf["pth2r_status"] == "NC"]
#print(list(subdf.index))

#print(subdf)
"""
#lets try again
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

#node initiation
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

#fill the leaves with fasta raw data
#returns the list of leaf node objects
leaflist = mf.get_leaf_list(n0)
for leaf in leaflist:
    #read the raw data. raw data directory example:
    #"/cta/users/ofkonar/work/resources/class_B1/canonical/calcr_orthologs_msa_pruned.fasta"
    rawdata = mf.read_fasta(fastadir + leaf.name + "_orthologs_msa_pruned.fasta")
    rawdatadf = mf.fasta_to_dataframe(rawdata)
    leaf.rawdata = rawdatadf

#fill the nodes with data, consensus of all leaves under the node
#also fill the nodes with metadata inner_check, where the positions are consensus
#get list of nodes
nodelist = mf.get_node_list(n0)
for node in nodelist:
    #for internal nodes, data is consensus of combination of leaves under the node
    if not mf.is_leaf(node):
        #get raw data of all leaves under
        noderawdatadict = mf.get_leaf_rawdata(node)
        #combine the sequences of leaves
        noderawdatalist = [noderawdatadict[key] for key in noderawdatadict]
        noderawdatadf = pd.concat(noderawdatalist, axis=1)
        #get the consensus of combined sequences
        nodeconsensus = mf.consensus(noderawdatadf, 0.9)
        #add node name to the data (to prevent confusion down the line)
        suffix = "_" + node.name
        nodeconsensus = nodeconsensus.add_suffix(suffix)
        #assign data to the node
        node.data = nodeconsensus
        node_con_psts = mf.df_equalto(nodeconsensus, 2, "C")
        node_metadata = {"internal_check" : node_con_psts}
        node.metadata = node_metadata

    else: #if the node is a leaf, we don't have to combine multiple leaves, rest is the same with inner nodes
        noderawdatadict = mf.get_leaf_rawdata(node)
        noderawdatadf = noderawdatadict[node.name]
        nodeconsensus = mf.consensus(noderawdatadf, 0.9)
        suffix = "_" + node.name
        nodeconsensus = nodeconsensus.add_suffix(suffix)
        node.data = nodeconsensus
        node_con_psts = mf.df_equalto(nodeconsensus, 2, "C")
        node_metadata = {"internal_check" : node_con_psts}
        node.metadata = node_metadata

#define ancestral trails
#ancestral trail of each receptor
glp2rlist = [n0, n1, glp2r]
glp1rlist = [n0, n1, n3, glp1r]
giprlist = [n0, n1, n3, n6, gipr]
gcgrlist = [n0, n1, n3, n6, gcgr]
sctrlist = [n0, n2, n4, sctr]
pacrlist = [n0, n2, n4, n7, pacr]
vipr2list = [n0, n2, n4, n7, n10, vipr2]
ghrhrlist = [n0, n2, n4, n7, n10, n13, ghrhr]
vipr1list = [n0, n2, n4, n7, n10, n13, vipr1]
pth1rlist = [n0, n2, n5, n8, pth1r]
pth2rlist = [n0, n2, n5, n8, pth2r]
crfr2list = [n0, n2, n5, n9, n11, crfr2]
crfr1list = [n0, n2, n5, n9, n11, crfr1]
calcrlist = [n0, n2, n5, n9, n12, calcr]
calrllist = [n0, n2, n5, n9, n12, calrl]

ancestrydict = {"glp2r" : glp2rlist, "glp1r" : glp1rlist, "gipr" : giprlist, "gcgr" : gcgrlist, "sctr" : sctrlist, "pacr" : pacrlist, "vipr2" : vipr2list, "ghrhr" : ghrhrlist, 
"vipr1" : vipr1list, "pth1r" : pth1rlist, "pth2r" : pth2rlist, "crfr2" : crfr2list, "crfr1" : crfr1list, "calcr" : calcrlist, "calrl" : calrllist}

for key in ancestrydict: #we will traverse from root to the node (follow the ancestral trail)
    ancestrylist = ancestrydict[key] #get the ancestral trail
    print("current ancestral trail leaf: ", key)
    seen = [] #nodes that are alraady seen while traversing the trail
    for i in range(len(ancestrylist)): #iterate over the trail
        current_node = ancestrylist[i] #current node that we want to compare with its ancestors
        print("current node to compare: ", current_node.name)
        if seen: #if there is any node in the seen list
            current_consensus = current_node.metadata["internal_check"] #get the consensus positions of the current node
            differences = [] #node vs ancestor list
            for j in seen: #iterate over seen nodes (ancestors)
                diff = mf.df_diff(current_node.data, 0, j.data, 0) #get the node vs ancestor difference
                checked_diff = [ele for ele in diff if ele in current_consensus] #get position that are different and consensus
                differences.append(checked_diff) #add them to the differences list
            total_ancestral_diff = list(set.intersection(*[set(x) for x in differences])) #get the intersections of all differences
            total_ancestral_diff.sort()
            current_node.metadata["total_ancestral_diff"] = total_ancestral_diff #assign the total ancestral difference to the metadata
            print(total_ancestral_diff)            
        seen.append(current_node)
    print("\n")

#compare each node with its sister
for node in nodelist:
    if not mf.is_leaf(node): #if node is a leaf it has no children to compare each other to
        left_children = node.left #get the left child
        left_consensus = left_children.metadata["internal_check"] #get the left consensus
        left_data = left_children.data #get the left data
        right_children = node.right
        right_consensus = right_children.metadata["internal_check"]
        right_data = right_children.data
        left_vs_right = mf.df_diff(left_data, 0, right_data, 0) #compare left and right
        leftdif = [ele for ele in left_vs_right if ele in left_consensus] #take positions that are different and consensus
        left_children.metadata["sister_difference"] = leftdif #add to the metadata
        right_vs_left = mf.df_diff(right_data, 0, left_data, 0)
        rightdif = [ele for ele in right_vs_left if ele in right_consensus]
        right_children.metadata["sister_difference"] = rightdif

for node in nodelist:
    if node != n0:
        print("node name: ", node.name)
        print([ele for ele in node.metadata["total_ancestral_diff"] if ele in node.metadata["sister_difference"]])
        print("\n")

print("My program took", time.time() - start_time, "seconds to run")