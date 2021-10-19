from functions import my_functions as mf
import pandas as pd
import time
from os import listdir
from os.path import isfile, join
start_time = time.time()

structurepath = "/cta/users/ofkonar/work/structures/selected/"
aaletterpath = pd.read_csv("/cta/users/ofkonar/work/resources/amino_acid_codes.csv")
structureinfo = pd.read_csv("/cta/users/ofkonar/work/resources/selected_structure.csv")
csvpath = "/cta/users/ofkonar/work/results/csvs/"
#rrcspath = "/cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py"
#for structure in files:
	#print("python " + rrcspath + " " + structurepath + structure)

"""
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/CALCR_6niy_Gs_Calcitonin_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/CRFR1_6pb0_Gs-i_UCN_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/CALRL_6e3y_Gs_RAMP1_CGRP_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/PTH1R_6fj3_PTH_Inactive.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/CRFR2_6pb1_Gs_UCN_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/CALRL_6uva_Gs_RAMP3_ADM2_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/VIPR1_6vn7_Gs_PACAP_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/CALRL_7knu_RAMP1_CGRP_Inactive.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/GLP1R_6ln2_NAM_Inactive.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/GHRHR_7cz5_Gs_Somatoliberin_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/GCGR_6lml_Gi_Glucagon_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/GLP1R_6x18_Gs_GLP1_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/CALRL_6uun_Gs_RAMP2_ADM_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/GCGR_5yqz_Partialagonist_Inactive.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/PTH1R_6nbf_Gs_PTHanalog_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/CRFR1_4k5y_T4lyzo_Antagonist_Inactive.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/CALRL_6uus_Gs_RAMP3_ADM_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/GCGR_6wpw_Gs_Glucagonder_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/GLP1R_7lci_Gs_Agonist_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/PACR_6m1i_Gs_PACAP_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/GLP2R_7d68_Gs_Glucagon_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/CRFR1_6p9x_Gs_Corticoliberin_Active.pdb
python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/selected/SCTR_6wzg_Gs_Secretin_Active.pdb
"""
"""
files = [f for f in listdir(structurepath) if isfile(join(structurepath, f))]
scorefiles = [x for x in files if "score" in x and "csv" if "csv" not in x]
for file in scorefiles:
	score = pd.read_csv(structurepath + file, delim_whitespace=True, header=None, engine="python")
	result = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"])
	dims = score.shape
	for i in range(dims[0]):
		line = score.iloc[i,]
		res1 = line[0]
		chain1 = res1.split(":")[0]
		aa1 = res1.split(":")[1]
		pos1 = aa1.split("_")[0]
		threecode1 = aa1.split("_")[1].lower()
		oneletter1 = aaletterpath.loc[aaletterpath["three_letter"] == threecode1, ["one_letter"]].values[0][0]
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

files = [f for f in listdir(structurepath) if isfile(join(structurepath, f))]
scorefiles = [x for x in files if "score.csv" in x]

for file in scorefiles:
	total = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"])
	prtname = file.split("_")[0]
	strname = file.split("_")[1].upper()
	targetcomparison = pd.read_csv(csvpath + prtname.lower() + "_binary_comparisons.csv")
	seqcolname = prtname.lower() + "_seq"
	fourteens = targetcomparison.loc[targetcomparison["count"] == 14, seqcolname].tolist()
	targetfile = pd.read_csv(structurepath + file)
	dims = targetfile.shape
	chain = structureinfo.loc[structureinfo["PDB_code"] == strname, "Chain_name"].values[0]
	for i in range(dims[0]):
		line = targetfile.iloc[i,]
		if line[0] == chain:
			if line[1] in fourteens:
				total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : line[4]}, ignore_index = True)
		elif line[2] == chain:
			if line[3] in fourteens:
				total = total.append({"chain_1" : line[0], "residue_1" : line[1], "chain_2" : line[2], "residue_2" : line[3], "RRCS" : line[4]}, ignore_index = True)
	total.to_csv("/cta/users/ofkonar/work/results/rrcs/selected/" + prtname + "_" + strname + "_" + "rrcs.csv", index = False)


print("My program took", time.time() - start_time, "seconds to run")