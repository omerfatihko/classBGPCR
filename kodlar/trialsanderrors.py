import pandas as pd
from os import listdir
from os.path import isfile, join
from collections import Counter
from functions import my_functions as mf

########################################################################################################################
#glp1r_str_path = "/cta/users/ofkonar/work/structures/GLP1R"
#glp1r_str_info = pd.read_csv("/cta/users/ofkonar/work/resources/glp1r_structure.csv", header = 0, index_col = 0)
#glp1r_files = [f for f in listdir(glp1r_str_path) if isfile(join(glp1r_str_path, f))]
#glp1r_scores = [i for i in glp1r_files if "score" in i]

#for file in glp1r_scores:
#	structure_name = file.split(".")[0].upper()
#	state = glp1r_str_info.loc[structure_name,'Activity']
#	a = pd.read_csv(glp1r_str_path + "/" + file, delim_whitespace=True, header=None, engine="python")
#	print(structure_name + "_" + state)
#	list_of_rows = []
#	for i in range(a.shape[0]):
#		row = a.iloc[i].tolist()
#		if any("362_THR" in ele for ele in row[0:2]):
#			list_of_rows.append(row)
#	df = pd.DataFrame(list_of_rows)
#	print(df)
########################################################################################################################
#get the column names
#colnames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sacc", "stitle"]
#read the blast results file
#ghrhrorthologs = pd.read_csv("/cta/users/ofkonar/work/results/secretin_subfamily/ghrhr/ghrhr_blast_results.csv", header = None, names = colnames)
#open the files where clade is kept and where we will write the ids and descriptions
#with open("/cta/users/ofkonar/work/tae/ghrhr_tree_weird_clade.txt", "r") as readfile, open("/cta/users/ofkonar/work/tae/ghrhr_clade_with_desc.txt", "w") as writefile:
#	#read the file line by line
#	for line in readfile:
#		#get the description of the protein
#		description = ghrhrorthologs.loc[ghrhrorthologs["sseqid"] == line[:-1], ["stitle"]].values[0][0]
#		#write the description and id to the file
#		writefile.write(description + " " + line)
########################################################################################################################
#colnames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sacc", "stitle"]
#read blast results
#ghrhrresults = pd.read_csv("/cta/users/ofkonar/work/results/secretin_subfamily/ghrhr/ghrhr_blast_results.csv", header = None, names = colnames)
#we get subject sequence ids (unique for every UniProt sequence, https://www.uniprot.org/help/fasta-headers)
#ghrhrids = ghrhrresults["sseqid"]
#get indices of human proteins
#indices = []
#for i, elem in enumerate(ghrhrids):
#    if 'HUMAN' in elem:
#        indices.append(i)
#target= indices[2]
#get all protein sequences until third human sequence (third human seq will act as outgroup)
#ghrhrselected = ghrhrids[0:target+1].tolist()
#get fastas from unified fasta as a list
#ghrhrlist = mf.filter_fasta(ghrhrselected, "/cta/users/ofkonar/work/resources/fasta/unified_fasta.fasta")
#write the fasta list to the path
#mf.write_fasta("/cta/users/ofkonar/work/tae/ghrhr_filtered_fastas_blast_third.fasta", ghrhrlist)
########################################################################################################################
#get the column names
colnames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sacc", "stitle"]
#read the blast results file
vipr1orthologs = pd.read_csv("/cta/users/ofkonar/work/results/secretin_subfamily/vipr1/vipr1_blast_results.csv", header = None, names = colnames)
#open the files where clade is kept and where we will write the ids and descriptions
with open("/cta/users/ofkonar/work/tae/vipr1_tree_weird_clade.txt", "r") as readfile, open("/cta/users/ofkonar/work/tae/vipr1_clade_with_desc.txt", "w") as writefile:
#	#read the file line by line
	for line in readfile:
		#get the description of the protein
		description = vipr1orthologs.loc[vipr1orthologs["sseqid"] == line[:-1], ["stitle"]].values[0][0]
		#write the description and id to the file
		writefile.write(description + " " + line)