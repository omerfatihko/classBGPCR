import pandas as pd
import time
start_time = time.time()

taxa_blast = pd.read_csv("/cta/users/ofkonar/work/results/vipr2_blast_results.csv", header = None, index_col = None)
proteins = taxa_blast[1].tolist()
#print(proteins)
orthologs = pd.DataFrame(columns=["protein_id", "header", "organism_id"])

with open("/cta/users/ofkonar/work/tae/taxa_tips.txt", "r") as file:
	for line in file:
		line = line.rstrip("\n")
		ind = proteins.index(line)
		row = taxa_blast.iloc[ind]
		header = row[13]
		parts = header.split(" ")
		imp = [x for x in parts if "OX=" in x][0]
		organism_id = imp.split("=")[1]
		orthologs = orthologs.append({"protein_id" : line, "header" : header, "organism_id" : organism_id}, ignore_index = True)
orthologs.to_csv("/cta/users/ofkonar/work/results/vipr2_orthologs.csv", index = False)
print("My program took", time.time() - start_time, "seconds to run")
