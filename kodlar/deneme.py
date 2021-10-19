import pandas as pd
from os import listdir
from os.path import isfile, join
from collections import Counter

gcgr_str_path = "/cta/users/ofkonar/work/structures/GCGR"
gcgr_str_info = pd.read_csv("/cta/users/ofkonar/work/resources/gcgr_structure.csv", header = 0, index_col = 0)
print(gcgr_str_info)

gcgr_files = [f for f in listdir(gcgr_str_path) if isfile(join(gcgr_str_path, f))]
gcgr_scores = [i for i in gcgr_files if "score" in i]
print(gcgr_scores)

for file in gcgr_scores:
	structure_name = file.split(".")[0].upper()
	print(structure_name)
	chain = gcgr_str_info.loc[structure_name,'Chain_name']
	print(chain)