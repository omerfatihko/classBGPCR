from functions import my_functions as mf
import pandas as pd
import time
from os import listdir
from os.path import isfile, join
start_time = time.time()
pd.set_option("display.max_rows", None, "display.max_columns", None) #  
"""READ ME!
These structures are selected from GPCRdb database at 06.09.2021"""
print("These structures are selected from GPCRdb database at 06.09.2021")
structurepath = "/cta/users/ofkonar/work/structures/all_structures/"
aaletterpath = pd.read_csv("/cta/users/ofkonar/work/resources/amino_acid_codes.csv")
structureinfo = pd.read_csv("/cta/users/ofkonar/work/resources/all_structures_info.csv")
csvpath = "/cta/users/ofkonar/work/results/csvs/"
rrcspath = "/cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py"

#files = [f for f in listdir(structurepath) if isfile(join(structurepath, f))] 

#print the commands to run RRCS code
#for structure in files:
#	print("python " + rrcspath + " " + structurepath + structure)

#RRCS commands
"""
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6vcb.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/7c2e.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6uus.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6x19.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6fj3.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/4l6r.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6orv.cif
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/7lci.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6xox.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/5vew.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6x18.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6vn7.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6uun.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/7cz5.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/5vex.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6e3y.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/7d3s.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6m1i.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6wzg.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/7lck.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6kk7.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/5xf1.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6nbh.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/5uz7.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/5xez.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/7knu.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6nbf.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6uva.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6lml.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6lpb.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6ln2.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6lmk.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6kjv.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6p9y.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/7d68.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6pb0.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6wpw.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/5ee7.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/7knt.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6b3j.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6nbi.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6kk1.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6wi9.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/4k5y.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/7lcj.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6m1h.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6x1a.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6pb1.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6p9x.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6niy.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/4z9g.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/6whc.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/all_structures/5yqz.pdb
"""
"""
#transform the results to a more practival form

files = [f for f in listdir(structurepath) if isfile(join(structurepath, f))]
scorefiles = [x for x in files if "score" in x and "csv" if "csv" not in x]

for file in scorefiles:
	score = pd.read_csv(structurepath + file, delim_whitespace=True, header=None, engine="python") #variable white space seperated files
	result = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"]) #open the dataframe to write transformed data
	dims = score.shape
	for i in range(dims[0]):
		line = score.iloc[i,]
		res1 = line[0]
		chain1 = res1.split(":")[0]
		aa1 = res1.split(":")[1]
		pos1 = aa1.split("_")[0]
		threecode1 = aa1.split("_")[1].lower()
		oneletter1 = aaletterpath.loc[aaletterpath["three_letter"] == threecode1, ["one_letter"]].values[0][0] #use aa letter code info 
		ele1 = oneletter1 + pos1
		res2 = line[1]
		chain2 = res2.split(":")[0]
		aa2 = res2.split(":")[1]
		pos2 = aa2.split("_")[0]
		threecode2 = aa2.split("_")[1].lower()
		oneletter2 = aaletterpath.loc[aaletterpath["three_letter"] == threecode2, ["one_letter"]].values[0][0]
		ele2 = oneletter2 + pos2
		rrcs = line[2]
		result = result.append({"chain_1" : chain1, "residue_1" : ele1, "chain_2" : chain2, "residue_2" : ele2, "RRCS" : rrcs}, ignore_index = True)
	result.to_csv(structurepath + file + ".csv", index = False)	
"""

#read the results for each GPCR, get global consensus and special residue interactions. Later, we will form a network with them and 
#supplement it with possible PP-interaction data.
#first we will work on active structures since most of the GPCRs have active structures while not all have inactive structures

files = [f for f in listdir(structurepath) if isfile(join(structurepath, f))]
scorefiles = [x for x in files if "score.csv" in x] #get transformed results

#calcr
"""
#get calcr active structures
calcrbinarycomp = pd.read_csv("/cta/users/ofkonar/work/results/csvs/calcr_binary_comparisons.csv") #read calcr binary comparison results
calcrspecres = calcrbinarycomp.loc[calcrbinarycomp["count"] == 14, "calcr_seq"].tolist() #get special residues with 14 counts
calcrglobalres = calcrbinarycomp.loc[calcrbinarycomp["global_status"] == "C", "calcr_seq"].tolist() #get global consensus residues for calcr
calcrselected = set(calcrspecres + calcrglobalres) #only interactions that involve either global or special residues will be selected
#calcractint = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "mean_RRCS"]) #open a file to keep active interactions
calcrpdbs = structureinfo.loc[(structureinfo["Receptor"]=="CALCR") & (structureinfo["Activity"]=="A"),"PDB_code"].tolist() #get active calcr structure ids
calcrfiles = [x for x in scorefiles if x.split(".")[0].upper() in calcrpdbs] #get their RRCS results

for structure in calcrfiles:
	total = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"])# open a dataframe to write selected interactions
	targetfile = pd.read_csv(structurepath + structure) #read the structure result
	dims = targetfile.shape
	strname = structure.split(".")[0].upper()
	chain = structureinfo.loc[structureinfo["PDB_code"] == strname, "Receptor_chain"].values[0] #get the chain name of the receptor
	nanobody = structureinfo.loc[structureinfo["PDB_code"] == strname, "Nanobody"].values[0] #get the nanobody chain name
	resolution = structureinfo.loc[structureinfo["PDB_code"] == strname, "Resolution"].values[0]
	for i in range(dims[0]):
		line = targetfile.iloc[i,]
		if line[0] == chain:
			if line[2] != nanobody:
				if line[1] in calcrselected:
					total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : (line[4]/resolution)}, ignore_index = True)
		elif line[2] == chain:
			if line[0] != nanobody:
				if line[3] in calcrselected:
					total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : (line[4]/resolution)}, ignore_index = True)
	total.to_csv("/cta/users/ofkonar/work/results/rrcs/all/" + "calcr_" + strname + "_" + "rrcs_resadjusted.csv", index = False)

#combine active calcr structures
#get the calcr files
rrcsselecteds = [f for f in listdir("/cta/users/ofkonar/work/results/rrcs/all/") if isfile(join("/cta/users/ofkonar/work/results/rrcs/all/", f))]
calcrrrcs = [x for x in rrcsselecteds if "calcr" in x]
calcrdf = pd.DataFrame() #open an empty dataframe, we will append the dataframes to it

for file in calcrrrcs:
	df = pd.read_csv("/cta/users/ofkonar/work/results/rrcs/all/" + file) #read the file
	calcrdf = pd.concat([calcrdf, df], axis = 0, ignore_index = True) #concatanate with main df
temp = calcrdf.groupby(["chain_1","residue_1","chain_2","residue_2"]).mean() #group the interactions
#residue-residue interactions always follow the same rule; regardless of chain name, residue with smaller index
#comes first. Thanks to that fact we can use groupby function of Pandas package.
#We don't have to check whether residues changed place (order is always the same)
temp.reset_index(inplace = True) #reset the indexes back to columns 
temp.to_csv("/cta/users/ofkonar/work/results/rrcs/all/" + "calr_rrcs_combined.csv", index = False)
"""

#calrl

#RAMP1, single structure # lets add its own consensus residues

#get calrl active structures
calrlbinarycomp = pd.read_csv("/cta/users/ofkonar/work/results/csvs/calrl_binary_comparisons.csv") #read calrl binary comparison results
calrlspecres = calrlbinarycomp.loc[calrlbinarycomp["count"] == 14, "calrl_seq"].tolist() #get special residues with 14 counts
calrlglobalres = calrlbinarycomp.loc[calrlbinarycomp["global_status"] == "C", "calrl_seq"].tolist() #get global consensus residues for calrl
calrlconsensus = calrlbinarycomp.loc[calrlbinarycomp["calrl_frequency"] == 1, "calrl_seq"].tolist() #get local consensus residues for calrl
calrlselected = set(calrlspecres + calrlglobalres + calrlconsensus) #only interactions that involve either global or special residues will be selected
print(calrlselected, len(calrlselected))
calrlpdbs = structureinfo.loc[(structureinfo["Receptor"]=="CALRL") & (structureinfo["Activity"]=="A") & (structureinfo["RAMP1"]=="RAMP1"),"PDB_code"].tolist() #get active calrl structures with RAMP1
print(calrlpdbs)
calrlfiles = [x for x in scorefiles if x.split(".")[0].upper() in calrlpdbs] #get their RRCS results

for structure in calrlfiles:
	total = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"])# open a dataframe to write selected interactions
	targetfile = pd.read_csv(structurepath + structure) #read the structure result
	dims = targetfile.shape
	strname = structure.split(".")[0].upper()
	chain = structureinfo.loc[structureinfo["PDB_code"] == strname, "Receptor_chain"].values[0] #get the chain name of the receptor
	print(chain)
	nanobody = structureinfo.loc[structureinfo["PDB_code"] == strname, "Nanobody"].values[0] #get the nanobody chain name
	resolution = structureinfo.loc[structureinfo["PDB_code"] == strname, "Resolution"].values[0] #get the resolution, will be used to adjust strength of interaction
	for i in range(dims[0]):
		line = targetfile.iloc[i,]
		if line[0] == chain:
			if line[2] != nanobody:
				if line[1] in calrlselected:
					total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : (line[4]/resolution)}, ignore_index = True)
		elif line[2] == chain:
			if line[0] != nanobody:
				if line[3] in calrlselected:
					#print(line)
					total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : (line[4]/resolution)}, ignore_index = True)
	total.to_csv("/cta/users/ofkonar/work/results/rrcs/all/" + "calrl_RAMP1_" + strname + "_" + "rrcs_resadjusted.csv", index = False)

#combine active calrl-RAMP1 structures
#get the calrl files
rrcsselecteds = [f for f in listdir("/cta/users/ofkonar/work/results/rrcs/all/") if isfile(join("/cta/users/ofkonar/work/results/rrcs/all/", f))]
calrlrrcs = [x for x in rrcsselecteds if "calrl_RAMP1_" in x if "combined" not in x]
print(calrlrrcs)
calrldf = pd.DataFrame() #open an empty dataframe, we will append the dataframes to it

for file in calrlrrcs:
	df = pd.read_csv("/cta/users/ofkonar/work/results/rrcs/all/" + file) #read the file
	calrldf = pd.concat([calrldf, df], axis = 0, ignore_index = True) #concatanate with main df
temp = calrldf.groupby(["chain_1","residue_1","chain_2","residue_2"]).mean() #group the interactions
#residue-residue interactions always follow the same rule; regardless of chain name, residue with smaller index
#comes first. Thanks to that fact we can use groupby function of Pandas package.
#We don't have to check whether residues changed place (order is always the same)
temp.reset_index(inplace = True) #reset the indexes back to columns
print(temp)
temp.to_csv("/cta/users/ofkonar/work/results/rrcs/all/" + "calrl_RAMP1_rrcs_combined.csv", index = False)


#RAMP2, single structure
"""
#get calrl active structures
calrlbinarycomp = pd.read_csv("/cta/users/ofkonar/work/results/csvs/calrl_binary_comparisons.csv") #read calrl binary comparison results
calrlspecres = calrlbinarycomp.loc[calrlbinarycomp["count"] == 14, "calrl_seq"].tolist() #get special residues with 14 counts
calrlglobalres = calrlbinarycomp.loc[calrlbinarycomp["global_status"] == "C", "calrl_seq"].tolist() #get global consensus residues for calrl
calrlselected = set(calrlspecres + calrlglobalres) #only interactions that involve either global or special residues will be selected
calrlpdbs = structureinfo.loc[(structureinfo["Receptor"]=="CALRL") & (structureinfo["Activity"]=="A") & (structureinfo["RAMP2"]=="RAMP2"),"PDB_code"].tolist() #get active calrl structures with RAMP2
print(calrlpdbs)
calrlfiles = [x for x in scorefiles if x.split(".")[0].upper() in calrlpdbs] #get their RRCS results

for structure in calrlfiles:
	total = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"])# open a dataframe to write selected interactions
	targetfile = pd.read_csv(structurepath + structure) #read the structure result
	dims = targetfile.shape
	strname = structure.split(".")[0].upper()
	chain = structureinfo.loc[structureinfo["PDB_code"] == strname, "Receptor_chain"].values[0] #get the chain name of the receptor
	nanobody = structureinfo.loc[structureinfo["PDB_code"] == strname, "Nanobody"].values[0] #get the nanobody chain name
	resolution = structureinfo.loc[structureinfo["PDB_code"] == strname, "Resolution"].values[0] #get the resolution, will be used to adjust strength of interaction
	for i in range(dims[0]):
		line = targetfile.iloc[i,]
		if line[0] == chain:
			if line[2] != nanobody:
				if line[1] in calrlselected:
					total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : (line[4]/resolution)}, ignore_index = True)
		elif line[2] == chain:
			if line[0] != nanobody:
				if line[3] in calrlselected:
					total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : (line[4]/resolution)}, ignore_index = True)
	total.to_csv("/cta/users/ofkonar/work/results/rrcs/all/" + "calrl_RAMP2_" + strname + "_" + "rrcs_resadjusted.csv", index = False)

#combine active calrl-RAMP2 structures
#get the calrl files
rrcsselecteds = [f for f in listdir("/cta/users/ofkonar/work/results/rrcs/all/") if isfile(join("/cta/users/ofkonar/work/results/rrcs/all/", f))]
calrlrrcs = [x for x in rrcsselecteds if "calrl_RAMP2_" in x]
calrldf = pd.DataFrame() #open an empty dataframe, we will append the dataframes to it

for file in calrlrrcs:
	df = pd.read_csv("/cta/users/ofkonar/work/results/rrcs/all/" + file) #read the file
	calrldf = pd.concat([calrldf, df], axis = 0, ignore_index = True) #concatanate with main df
temp = calrldf.groupby(["chain_1","residue_1","chain_2","residue_2"]).mean() #group the interactions
#residue-residue interactions always follow the same rule; regardless of chain name, residue with smaller index
#comes first. Thanks to that fact we can use groupby function of Pandas package.
#We don't have to check whether residues changed place (order is always the same)
temp.reset_index(inplace = True) #reset the indexes back to columns 
temp.to_csv("/cta/users/ofkonar/work/results/rrcs/all/" + "carl_RAMP2_rrcs_combined.csv", index = False)
"""

#RAMP3
"""
#get calrl active structures
calrlbinarycomp = pd.read_csv("/cta/users/ofkonar/work/results/csvs/calrl_binary_comparisons.csv") #read calrl binary comparison results
calrlspecres = calrlbinarycomp.loc[calrlbinarycomp["count"] == 14, "calrl_seq"].tolist() #get special residues with 14 counts
calrlglobalres = calrlbinarycomp.loc[calrlbinarycomp["global_status"] == "C", "calrl_seq"].tolist() #get global consensus residues for calrl
calrlselected = set(calrlspecres + calrlglobalres) #only interactions that involve either global or special residues will be selected
print(len(calrlselected))
calrlpdbs = structureinfo.loc[(structureinfo["Receptor"]=="CALRL") & (structureinfo["Activity"]=="A") & (structureinfo["RAMP3"]=="RAMP3"),"PDB_code"].tolist() #get active calrl structures with RAMP3
print(calrlpdbs)
calrlfiles = [x for x in scorefiles if x.split(".")[0].upper() in calrlpdbs] #get their RRCS results

for structure in calrlfiles:
	total = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"])# open a dataframe to write selected interactions
	targetfile = pd.read_csv(structurepath + structure) #read the structure result
	dims = targetfile.shape
	strname = structure.split(".")[0].upper()
	chain = structureinfo.loc[structureinfo["PDB_code"] == strname, "Receptor_chain"].values[0] #get the chain name of the receptor
	nanobody = structureinfo.loc[structureinfo["PDB_code"] == strname, "Nanobody"].values[0] #get the nanobody chain name
	resolution = structureinfo.loc[structureinfo["PDB_code"] == strname, "Resolution"].values[0] #get the resolution, will be used to adjust strength of interaction
	for i in range(dims[0]):
		line = targetfile.iloc[i,]
		if line[0] == chain:
			if line[2] != nanobody:
				if line[1] in calrlselected:
					total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : (line[4]/resolution)}, ignore_index = True)
		elif line[2] == chain:
			if line[0] != nanobody:
				if line[3] in calrlselected:
					total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : (line[4]/resolution)}, ignore_index = True)
	total.to_csv("/cta/users/ofkonar/work/results/rrcs/all/" + "calrl_RAMP3_" + strname + "_" + "rrcs_resadjusted.csv", index = False)

#combine active calrl-RAMP3 structures
#get the calrl files
rrcsselecteds = [f for f in listdir("/cta/users/ofkonar/work/results/rrcs/all/") if isfile(join("/cta/users/ofkonar/work/results/rrcs/all/", f))]
calrlrrcs = [x for x in rrcsselecteds if "calrl_RAMP3_" in x]
calrldf = pd.DataFrame() #open an empty dataframe, we will append the dataframes to it

for file in calrlrrcs:
	df = pd.read_csv("/cta/users/ofkonar/work/results/rrcs/all/" + file) #read the file
	calrldf = pd.concat([calrldf, df], axis = 0, ignore_index = True) #concatanate with main df
temp = calrldf.groupby(["chain_1","residue_1","chain_2","residue_2"]).mean() #group the interactions
#residue-residue interactions always follow the same rule; regardless of chain name, residue with smaller index
#comes first. Thanks to that fact we can use groupby function of Pandas package.
#We don't have to check whether residues changed place (order is always the same)
temp.reset_index(inplace = True) #reset the indexes back to columns 
temp.to_csv("/cta/users/ofkonar/work/results/rrcs/all/" + "carl_RAMP3_rrcs_combined.csv", index = False)
"""

#crfr1
"""
#get crfr1 active structures
crfr1binarycomp = pd.read_csv("/cta/users/ofkonar/work/results/csvs/crfr1_binary_comparisons.csv") #read crfr1 binary comparison results
crfr1specres = crfr1binarycomp.loc[crfr1binarycomp["count"] == 14, "crfr1_seq"].tolist() #get special residues with 14 counts
crfr1globalres = crfr1binarycomp.loc[crfr1binarycomp["global_status"] == "C", "crfr1_seq"].tolist() #get global consensus residues for crfr1
crfr1selected = set(crfr1specres + crfr1globalres) #only interactions that involve either global or special residues will be selected
print(len(crfr1selected))
crfr1pdbs = structureinfo.loc[(structureinfo["Receptor"]=="CRFR1") & (structureinfo["Activity"]=="A"),"PDB_code"].tolist() #get active crfr1 structures
print(crfr1pdbs)
crfr1files = [x for x in scorefiles if x.split(".")[0].upper() in crfr1pdbs] #get their RRCS results

for structure in crfr1files:
	total = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"])# open a dataframe to write selected interactions
	targetfile = pd.read_csv(structurepath + structure) #read the structure result
	dims = targetfile.shape
	strname = structure.split(".")[0].upper()
	chain = structureinfo.loc[structureinfo["PDB_code"] == strname, "Receptor_chain"].values[0] #get the chain name of the receptor
	nanobody = structureinfo.loc[structureinfo["PDB_code"] == strname, "Nanobody"].values[0] #get the nanobody chain name
	resolution = structureinfo.loc[structureinfo["PDB_code"] == strname, "Resolution"].values[0] #get the resolution, will be used to adjust strength of interaction
	for i in range(dims[0]):
		line = targetfile.iloc[i,]
		if line[0] == chain:
			if line[2] != nanobody:
				if line[1] in crfr1selected:
					total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : (line[4]/resolution)}, ignore_index = True)
		elif line[2] == chain:
			if line[0] != nanobody:
				if line[3] in crfr1selected:
					total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : (line[4]/resolution)}, ignore_index = True)
	total.to_csv("/cta/users/ofkonar/work/results/rrcs/all/" + "crfr1_" + strname + "_" + "rrcs_resadjusted.csv", index = False)

#combine active crfr1 structures
#get the crfr1 files
rrcsselecteds = [f for f in listdir("/cta/users/ofkonar/work/results/rrcs/all/") if isfile(join("/cta/users/ofkonar/work/results/rrcs/all/", f))]
crfr1rrcs = [x for x in rrcsselecteds if "crfr1_" in x]
crfr1df = pd.DataFrame() #open an empty dataframe, we will append the dataframes to it

for file in crfr1rrcs:
	df = pd.read_csv("/cta/users/ofkonar/work/results/rrcs/all/" + file) #read the file
	crfr1df = pd.concat([crfr1df, df], axis = 0, ignore_index = True) #concatanate with main df
temp = crfr1df.groupby(["chain_1","residue_1","chain_2","residue_2"]).mean() #group the interactions
#residue-residue interactions always follow the same rule; regardless of chain name, residue with smaller index
#comes first. Thanks to that fact we can use groupby function of Pandas package.
#We don't have to check whether residues changed place (order is always the same)
temp.reset_index(inplace = True) #reset the indexes back to columns 
temp.to_csv("/cta/users/ofkonar/work/results/rrcs/all/" + "crfr1_rrcs_combined.csv", index = False)
"""

#crfr2 ile devam
print("My program took", time.time() - start_time, "seconds to run")