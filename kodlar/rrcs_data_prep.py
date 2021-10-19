import pandas as pd
from os import listdir
from os.path import isfile, join
from collections import Counter

#rrcs_path = "/cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py"
gcgr_str_path = "/cta/users/ofkonar/work/structures/GCGR"
gcgr_str_info = pd.read_csv("/cta/users/ofkonar/work/resources/gcgr_structure.csv", header = 0, index_col = 0)
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GCGR/4l6r.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GCGR/5xez.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GCGR/5xf1.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GCGR/6lmk.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GCGR/5ee7.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GCGR/6wpw.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GCGR/5yqz.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GCGR/6whc.pdb
glp1r_str_path = "/cta/users/ofkonar/work/structures/GLP1R"
glp1r_str_info = pd.read_csv("/cta/users/ofkonar/work/resources/glp1r_structure.csv", header = 0, index_col = 0)
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//5vex.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//6x18.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//6orv.cif
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//6xox.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//7lci.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//5vew.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//6vcb.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//6x19.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//6x1a.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//6b3j.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//7lcj.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//6kk1.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//6ln2.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//6kjv.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//7lck.pdb
#python /cta/users/ofkonar/work/kodlar/RRCS-master/RRCS.py /cta/users/ofkonar/work/structures/GLP1R//6kk7.pdb

#for file in glp1r_files:
#	print("python " + rrcs_path + " " + glp1r_str_path + "/" + file)

gcgr_files = [f for f in listdir(gcgr_str_path) if isfile(join(gcgr_str_path, f))]
gcgr_scores = [i for i in gcgr_files if "score" in i]

glp1r_files = [f for f in listdir(glp1r_str_path) if isfile(join(glp1r_str_path, f))]
glp1r_scores = [i for i in glp1r_files if "score" in i]

gcgr_scores_active = {}
gcgr_scores_inactive = {}
gcgr_length = 477
gcgr_active_structures = []
gcgr_inactive_structures = []

for file in gcgr_scores:
	structure_name = file.split(".")[0].upper()
	a = pd.read_csv(gcgr_str_path + "/" + file, delim_whitespace=True, header=None, engine="python")
	dims = a.shape
	b = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"])
	chain = gcgr_str_info.loc[structure_name,'Chain_name']
	state = gcgr_str_info.loc[structure_name,'Activity']
	if state == "A":
		gcgr_active_structures.append(structure_name)
		for i in range(dims[0]):
			row = a.iloc[i]
			chain_1 = row[0].split(":")[0]
			chain_2 = row[1].split(":")[0]
			res_1 = row[0].split(":")[1]
			res_2 = row[1].split(":")[1]
			pos_1 = int(''.join(x for x in res_1 if x.isdigit()))
			pos_2 = int(''.join(x for x in res_2 if x.isdigit()))
			if chain_1 == chain and chain_2 == chain and pos_1 < (gcgr_length + 1) and pos_2 < (gcgr_length + 1):
				b = b.append({"chain_1" : chain_1, "residue_1" : res_1, "chain_2" : chain_2, "residue_2" : res_2, "RRCS" : row[2]}, ignore_index = True)
		gcgr_scores_active[structure_name] = b
	elif state == "I":
		gcgr_inactive_structures.append(structure_name)
		for i in range(dims[0]):
			row = a.iloc[i]
			chain_1 = row[0].split(":")[0]
			chain_2 = row[1].split(":")[0]
			res_1 = row[0].split(":")[1]
			res_2 = row[1].split(":")[1]
			pos_1 = int(''.join(x for x in res_1 if x.isdigit()))
			pos_2 = int(''.join(x for x in res_2 if x.isdigit()))
			if chain_1 == chain and chain_2 == chain and pos_1 < (gcgr_length + 1) and pos_2 < (gcgr_length + 1):
				b = b.append({"chain_1" : chain_1, "residue_1" : res_1, "chain_2" : chain_2, "residue_2" : res_2, "RRCS" : row[2]}, ignore_index = True)
		gcgr_scores_inactive[structure_name] = b

glp1r_scores_active = {}
glp1r_scores_inactive = {}
glp1r_length = 463
glp1r_active_structures = []
glp1r_inactive_structures = []

for file in glp1r_scores:
	structure_name = file.split(".")[0].upper()
	a = pd.read_csv(glp1r_str_path + "/" + file, delim_whitespace=True, header=None, engine="python")
	dims = a.shape
	b = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"])
	chain = glp1r_str_info.loc[structure_name,'Chain_name']
	state = glp1r_str_info.loc[structure_name,'Activity']
	if state == "A":
		glp1r_active_structures.append(structure_name)
		for i in range(dims[0]):
			row = a.iloc[i]
			chain_1 = row[0].split(":")[0]
			chain_2 = row[1].split(":")[0]
			res_1 = row[0].split(":")[1]
			res_2 = row[1].split(":")[1]
			pos_1 = int(''.join(x for x in res_1 if x.isdigit()))
			pos_2 = int(''.join(x for x in res_2 if x.isdigit()))
			if chain_1 == chain and chain_2 == chain and pos_1 < (glp1r_length + 1) and pos_2 < (glp1r_length + 1):
				b = b.append({"chain_1" : chain_1, "residue_1" : res_1, "chain_2" : chain_2, "residue_2" : res_2, "RRCS" : row[2]}, ignore_index = True)
		glp1r_scores_active[structure_name] = b
	elif state == "I":
		glp1r_inactive_structures.append(structure_name)
		for i in range(dims[0]):
			row = a.iloc[i]
			chain_1 = row[0].split(":")[0]
			chain_2 = row[1].split(":")[0]
			res_1 = row[0].split(":")[1]
			res_2 = row[1].split(":")[1]
			pos_1 = int(''.join(x for x in res_1 if x.isdigit()))
			pos_2 = int(''.join(x for x in res_2 if x.isdigit()))			
			if chain_1 == chain and chain_2 == chain and pos_1 < (glp1r_length + 1) and pos_2 < (glp1r_length + 1):
				b = b.append({"chain_1" : chain_1, "residue_1" : res_1, "chain_2" : chain_2, "residue_2" : res_2, "RRCS" : row[2]}, ignore_index = True)
		glp1r_scores_inactive[structure_name] = b

gcgr_all_pairs = []

for i in gcgr_scores_active:
	df = gcgr_scores_active[i]
	res_pairs = df[["residue_1", "residue_2"]].to_records(index = False).tolist()
	for pair in res_pairs:
		if pair not in gcgr_all_pairs:
			gcgr_all_pairs.append(pair)
for i in gcgr_scores_inactive:
	df = gcgr_scores_inactive[i]
	res_pairs = df[["residue_1", "residue_2"]].to_records(index = False).tolist()
	for pair in res_pairs:
		if pair not in gcgr_all_pairs:
			gcgr_all_pairs.append(pair)

glp1r_all_pairs = []

for i in glp1r_scores_active:
	df = glp1r_scores_active[i]
	res_pairs = df[["residue_1", "residue_2"]].to_records(index = False).tolist()
	for pair in res_pairs:
		if pair not in glp1r_all_pairs:
			glp1r_all_pairs.append(pair)
for i in glp1r_scores_inactive:
	df = glp1r_scores_inactive[i]
	res_pairs = df[["residue_1", "residue_2"]].to_records(index = False).tolist()
	for pair in res_pairs:
		if pair not in glp1r_all_pairs:
			glp1r_all_pairs.append(pair)

gcgr_active_pairs = pd.DataFrame(columns = gcgr_active_structures, index = gcgr_all_pairs)
gcgr_inactive_pairs = pd.DataFrame(columns = gcgr_inactive_structures, index = gcgr_all_pairs)

for column in gcgr_active_pairs:
	df = gcgr_scores_active[column]
	for ele in gcgr_all_pairs:
		if ((df['residue_1'] == ele[0]) & (df['residue_2'] == ele[1])).any():
			RRCS = df[(df['residue_1'] == ele[0]) & (df['residue_2'] == ele[1])]["RRCS"].values[0]
		else:
			RRCS = 0
		gcgr_active_pairs.at[ele,column] = RRCS

for column in gcgr_inactive_pairs:
	df = gcgr_scores_inactive[column]
	for ele in gcgr_all_pairs:
		if ((df['residue_1'] == ele[0]) & (df['residue_2'] == ele[1])).any():
			RRCS = df[(df['residue_1'] == ele[0]) & (df['residue_2'] == ele[1])]["RRCS"].values[0]
		else:
			RRCS = 0
		gcgr_inactive_pairs.at[ele,column] = RRCS

glp1r_active_pairs = pd.DataFrame(columns = glp1r_active_structures, index = glp1r_all_pairs)
glp1r_inactive_pairs = pd.DataFrame(columns = glp1r_inactive_structures, index = glp1r_all_pairs)

for column in glp1r_active_pairs:
	df = glp1r_scores_active[column]
	for ele in glp1r_all_pairs:
		if ((df['residue_1'] == ele[0]) & (df['residue_2'] == ele[1])).any():
			RRCS = df[(df['residue_1'] == ele[0]) & (df['residue_2'] == ele[1])]["RRCS"].values[0]
		else:
			RRCS = 0
		glp1r_active_pairs.at[ele,column] = RRCS

for column in glp1r_inactive_pairs:
	df = glp1r_scores_inactive[column]
	for ele in glp1r_all_pairs:
		if ((df['residue_1'] == ele[0]) & (df['residue_2'] == ele[1])).any():
			RRCS = df[(df['residue_1'] == ele[0]) & (df['residue_2'] == ele[1])]["RRCS"].values[0]
		else:
			RRCS = 0
		glp1r_inactive_pairs.at[ele,column] = RRCS

gcgr_active_pairs.to_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_RRCS_active.csv")
gcgr_inactive_pairs.to_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_RRCS_inactive.csv")
glp1r_active_pairs.to_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_RRCS_active.csv")
glp1r_inactive_pairs.to_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_RRCS_inactive.csv")