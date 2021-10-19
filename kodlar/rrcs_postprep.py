import pandas as pd
from os import listdir
from os.path import isfile, join
from collections import Counter
from functions import my_functions as mf

gcgr_active = pd.read_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_RRCS_active.csv", header = 0, index_col = 0)
gcgr_inactive = pd.read_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_RRCS_inactive.csv", header = 0, index_col = 0)
glp1r_active = pd.read_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_RRCS_active.csv", header = 0, index_col = 0)
glp1r_inactive = pd.read_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_RRCS_inactive.csv", header = 0, index_col = 0)

wootten_table = pd.read_csv("/cta/users/ofkonar/work/resources/residue_table.csv")
gcgr_column = wootten_table["GCGR"]
glp1r_column = wootten_table["GLP1R"]
residue_table = pd.read_csv("/cta/users/ofkonar/work/resources/residue_table.csv")
amino_acid_codes = pd.read_csv("/cta/users/ofkonar/work/resources/amino_acid_codes.csv")

gcgr_active_wootten = pd.DataFrame(columns = list(gcgr_active.columns))

for ind in list(gcgr_active.index):
	parts = ind.split(",")
	part1 = "".join(c for c in parts[0] if c not in "'()")
	part2 = "".join(c for c in parts[1] if c not in "'()")
	index1 = part1.split("_")[0]
	aa1 = part1.split("_")[1].lower()
	index2 = part2.split("_")[0].replace(" ","")
	aa2 = part2.split("_")[1].lower().replace(" ","")
	res1 = amino_acid_codes.loc[amino_acid_codes["three_letter"] == aa1].one_letter.values[0] + index1
	res2 = amino_acid_codes.loc[amino_acid_codes["three_letter"] == aa2].one_letter.values[0] + index2
	if res1 in set(gcgr_column):
		if res2 in set(gcgr_column):
			wootten1 = str(wootten_table.loc[wootten_table["GCGR"] == res1, "Wootten"].values[0])
			wootten2 = str(wootten_table.loc[wootten_table["GCGR"] == res2, "Wootten"].values[0])
			new_index = wootten1 + "," + wootten2
			gcgr_active_wootten.loc[new_index] = gcgr_active.loc[ind]

glp1r_active_wootten = pd.DataFrame(columns = list(glp1r_active.columns))

for ind in list(glp1r_active.index):
	parts = ind.split(",")
	part1 = "".join(c for c in parts[0] if c not in "'()")
	part2 = "".join(c for c in parts[1] if c not in "'()")
	index1 = part1.split("_")[0]
	aa1 = part1.split("_")[1].lower()
	index2 = part2.split("_")[0].replace(" ","")
	aa2 = part2.split("_")[1].lower().replace(" ","")
	res1 = amino_acid_codes.loc[amino_acid_codes["three_letter"] == aa1].one_letter.values[0] + index1
	res2 = amino_acid_codes.loc[amino_acid_codes["three_letter"] == aa2].one_letter.values[0] + index2
	if res1 in set(glp1r_column):
		if res2 in set(glp1r_column):
			wootten1 = str(wootten_table.loc[wootten_table["GLP1R"] == res1, "Wootten"].values[0])
			wootten2 = str(wootten_table.loc[wootten_table["GLP1R"] == res2, "Wootten"].values[0])
			new_index = wootten1 + "," + wootten2
			glp1r_active_wootten.loc[new_index] = glp1r_active.loc[ind]

shared_rows = [val for val in list(glp1r_active_wootten.index) if val in list(gcgr_active_wootten.index)]
gcgr_glp1r_combined_active_wootten = pd.concat([gcgr_active_wootten.loc[shared_rows], glp1r_active_wootten.loc[shared_rows]], axis = 1)
gcgr_glp1r_combined_active_wootten.to_csv("/cta/users/ofkonar/work/results/gcgr_glp1r_combined_active_wootten.csv")

gcgr_inactive_wootten = pd.DataFrame(columns = list(gcgr_inactive.columns))

for ind in list(gcgr_inactive.index):
	parts = ind.split(",")
	part1 = "".join(c for c in parts[0] if c not in "'()")
	part2 = "".join(c for c in parts[1] if c not in "'()")
	index1 = part1.split("_")[0]
	aa1 = part1.split("_")[1].lower()
	index2 = part2.split("_")[0].replace(" ","")
	aa2 = part2.split("_")[1].lower().replace(" ","")
	res1 = amino_acid_codes.loc[amino_acid_codes["three_letter"] == aa1].one_letter.values[0] + index1
	res2 = amino_acid_codes.loc[amino_acid_codes["three_letter"] == aa2].one_letter.values[0] + index2
	if res1 in set(gcgr_column):
		if res2 in set(gcgr_column):
			wootten1 = str(wootten_table.loc[wootten_table["GCGR"] == res1, "Wootten"].values[0])
			wootten2 = str(wootten_table.loc[wootten_table["GCGR"] == res2, "Wootten"].values[0])
			new_index = wootten1 + "," + wootten2
			gcgr_inactive_wootten.loc[new_index] = gcgr_inactive.loc[ind]

glp1r_inactive_wootten = pd.DataFrame(columns = list(glp1r_inactive.columns))

for ind in list(glp1r_inactive.index):
	parts = ind.split(",")
	part1 = "".join(c for c in parts[0] if c not in "'()")
	part2 = "".join(c for c in parts[1] if c not in "'()")
	index1 = part1.split("_")[0]
	aa1 = part1.split("_")[1].lower()
	index2 = part2.split("_")[0].replace(" ","")
	aa2 = part2.split("_")[1].lower().replace(" ","")
	res1 = amino_acid_codes.loc[amino_acid_codes["three_letter"] == aa1].one_letter.values[0] + index1
	res2 = amino_acid_codes.loc[amino_acid_codes["three_letter"] == aa2].one_letter.values[0] + index2
	if res1 in set(glp1r_column):
		if res2 in set(glp1r_column):
			wootten1 = str(wootten_table.loc[wootten_table["GLP1R"] == res1, "Wootten"].values[0])
			wootten2 = str(wootten_table.loc[wootten_table["GLP1R"] == res2, "Wootten"].values[0])
			new_index = wootten1 + "," + wootten2
			glp1r_inactive_wootten.loc[new_index] = glp1r_inactive.loc[ind]

shared_rows = [val for val in list(glp1r_inactive_wootten.index) if val in list(gcgr_inactive_wootten.index)]
gcgr_glp1r_combined_inactive_wootten = pd.concat([gcgr_inactive_wootten.loc[shared_rows], glp1r_inactive_wootten.loc[shared_rows]], axis = 1)
gcgr_glp1r_combined_inactive_wootten.to_csv("/cta/users/ofkonar/work/results/gcgr_glp1r_combined_inactive_wootten.csv")