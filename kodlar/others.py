from functions import my_functions as mf
import pandas as pd
from os import listdir
from os.path import isfile, join
from collections import Counter
import re
import time
start_time = time.time()
unified_fasta_path = "/cta/users/ofkonar/work/database/canonical_unified_fasta/unified_fasta.fasta"
others_fasta_path = "/cta/users/ofkonar/work/resources/class_B1/"
general_fasta_path = "/cta/users/ofkonar/work/resources/fasta/"
pd.set_option("display.max_columns", None) # "display.max_rows", None, 
prt_list = ["calcr", "calrl", "crfr1", "crfr2", "gcgr", "ghrhr", "gipr", "glp1r", "glp2r", "pacr",  "pth1r", "pth2r", "sctr", "vipr1", "vipr2"]
csv_path = "/cta/users/ofkonar/work/results/"
#
############################################################################################################################

#filter blast results
#for prt in prt_list:
	#set current proteins csv file path
	#currentpath = csv_path + prt

#	#get csv file containing blast results from the said path
#	files = [f for f in listdir(csv_path) if isfile(join(csv_path, f))]
#	target = [x for x in files if prt in x]
#	csvfile = csv_path + target[0]
	
	#read the blast results
#	selected = mf.select_csv(csvfile, 1)
	#filter the results and write them down
#	blastfastas = mf.filter_fasta(selected, unified_fasta_path)
#	mf.write_fasta(others_fasta_path + prt + "_filtered_fastas_blast.fasta", blastfastas)

#selected = mf.select_csv("/cta/users/ofkonar/work/results/glp1r_blast_results.csv", 1)
#for i in selected:
#	print(i)
#print("******************************************************************************************")
#print(len(selected))
#print("******************************************************************************************")
#blastfastas = mf.filter_fasta(selected, unified_fasta_path)

#print(len(blastfastas))
#print("******************************************************************************************")
#for i in blastfastas:
#	print(i)
#print("******************************************************************************************")
#mf.write_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/glp1r_filtered_fastas_blast.fasta", blastfastas)
############################################################################################################################

#filter all class B1 orthologs
#path of orthologs file
#orthologsfile = "/cta/users/ofkonar/work/results/csvs/class_B1_all_orthologs.csv"
#select orthologs
#selected = mf.select_csv(orthologsfile, 0)
#filter the results and write them down
#allorthologfastas = mf.filter_fasta(selected, unified_fasta_path)
#mf.write_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/class_B1_all_orthologs.fasta",allorthologfastas)
############################################################################################################################

#prune the gaps in class B1 orthologs
#read the fasta file
#classb1orths = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/class_B1_all_orthologs_msa.fasta")
#pruning step
#prunedclassb1 = mf.fasta_pruner(classb1orths)
#write down the file 
#mf.write_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/class_B1_all_orthologs_msa_pruned.fasta", prunedclassb1)
############################################################################################################################

#seperate class B1 proteins
currentcsvpath = "/cta/users/ofkonar/work/results/csvs/"
currentfastapath = "/cta/users/ofkonar/work/resources/class_B1/canonical/"
#get the csv files with seperated orthologs
#filename = [f for f in listdir(currentcsvpath) if isfile(join(currentcsvpath, f))]
#memberfiles = []
#for f in filename:
#	if ("_orthologs" in f) and ("all" not in f):
#		memberfiles.append(f)
#		print(f)

#for f in memberfiles:
	#select the orthologs
#	ortfile = currentcsvpath + f
#	print(ortfile)
#	selected = mf.select_csv(ortfile, 0)
	#filter them from the fasta list
#	fastas = mf.filter_fasta(selected, "/cta/users/ofkonar/work/resources/class_B1/canonical/class_B1_all_orthologs_msa_pruned.fasta")
	#write them down
#	write_path = currentfastapath + f.split("_")[0] + "_orthologs_msa_pruned.fasta"
#	mf.write_fasta(write_path, fastas)
############################################################################################################################

#get consensus for individual proteins (prune the gaps in them)
#get the fasta files
#filename = [f for f in listdir(currentfastapath) if isfile(join(currentfastapath, f))]
#targetfiles = [f for f in filename if "_pruned" in f if "B1" not in f]
#for target in targetfiles:
	#read the fasta file
#	temp = mf.read_fasta(currentfastapath + target)
	#prune the gaps
#	pruned = mf.fasta_pruner(temp)
	#turn pruned fasta to df
#	pruneddf = mf.fasta_to_dataframe(pruned)
	#get the consensus
#	consensus = mf.consensus(pruneddf, 0.9)
	#write down the consensus
#	consensus.to_csv(currentcsvpath + target.split("_")[0] + "_consensus.csv")
############################################################################################################################

#all comparisons

#read the Wootten file for later use
Wootten_table = pd.read_csv("/cta/users/ofkonar/work/resources/residue_table.csv")
#read the domain file for later use
domain_table = pd.read_csv("/cta/users/ofkonar/work/resources/class_B1_domains.csv", index_col = 0)

#get global consensus
#read the pruned global alignment
temp = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/class_B1_all_orthologs_msa_pruned.fasta")
#turn it into dataframe
tempdf = mf.fasta_to_dataframe(temp)
#get the column names (for later use)
seqnames = list(tempdf.columns)
#get the human sequence names
humanseqs = [x for x in seqnames if "HUMAN" in x] # if prt in x.lower()
#get uniprot ids for human proteins, keep them in a dictionary
prttoid = {}
for item in humanseqs:
	prtidentifier = item.split(" ")[0]
	prtid = prtidentifier.split("|")[1]
	prtname = prtidentifier.split("|")[2].split("_")[0]
	if prtname == "GLR":
		prtname = "GCGR"
	prttoid[prtname] = prtid

#get the consensus
globalconsensus = mf.consensus(tempdf, 0.9)
#write global consensus down
#globalconsensus.to_csv(currentcsvpath + "class_B1_global_consensus.csv") 
#add global_ prefix to the columns
globalconsensus = globalconsensus.add_prefix("global_") ###built upon this


#get the fasta files
filename = [f for f in listdir(currentfastapath) if isfile(join(currentfastapath, f))]
targetfiles = [f for f in filename if "_pruned" in f if "B1" not in f]

#open a dictionary to hold consensus dataframes
consensusdict = {}
#open a list to hold protein names, these names will be used as index for the dictionary
prtnames = []
for target in targetfiles:
	#get consensus of the target file
	#read the target fasta file
	targetfasta = mf.read_fasta(currentfastapath + target)
	#turn it into dataframe
	targetdf = mf.fasta_to_dataframe(targetfasta)
	#get the consensus
	targetconsensus = mf.consensus(targetdf, 0.9)
	#get prt name 
	targetname = target.split("_")[0]
	#add the prt name to the list
	prtnames.append(targetname)
	#add consensus df to the dictionary
	consensusdict[targetname] = targetconsensus

#compare the proteins
for prt in prtnames:
	print(prt)
	print("**************************************")	
	#get rest of the proteins to compare
	rest = list(prtnames)
	rest.remove(prt)
	#concatenate global consensus and local consensus, we will build upon this
	localconsensus = consensusdict[prt]
	localconsensus =localconsensus.add_prefix(prt + "_")
	result = pd.concat([globalconsensus, localconsensus], axis = 1)

	#get the human sequence of targeted protein, add it to the dataframe
	#get the uniprot id of target protein
	kilit = prttoid[prt.upper()]
	#get the sequence column with required id
	tarcolname = [x for x in seqnames if kilit in x] 
	tarcol = tempdf[tarcolname]
	#it is a one column dataframe so we squeeze it to a series
	tarcol = tarcol.squeeze()
	counter = 1
	Wootten = []
	domains = []
	#add residue index to residue name, create Wootten column
	#get the required protein column
	my_column = Wootten_table[prt.upper()]
	for items in tarcol.iteritems():
		if items[1] != "-":
			newres = items[1] + str(counter)
			tarcol.set_value(items[0], newres)
			#get domain info of the residue
			a = domain_table.loc[domain_table[prt.upper()] >=counter, prt.upper()]
			#print(counter)
			#print(a)
			domain = a.index.values[0]
			#print(domain)
			domains.append(domain)
			counter += 1
			if newres in set(my_column):
				#get the index of Wootten
				windex = my_column[my_column == newres].index[0]
				#append to the Wootten column
				Wootten.append(str(Wootten_table.iloc[windex]['Wootten']))
			else:
				#if res lacks Wootten numbering
				Wootten.append('-')
		else:
			#for gaps
			Wootten.append('-')
			domains.append("-")
	result.insert(result.shape[1], "Wootten", Wootten, True)
	result.insert(result.shape[1], "Domain", domains, True)
	result.insert(3, prt + "_seq", tarcol)
	
	dfsize = result.shape[1]
	#get blosum scores for each comparison
	for i in rest:
		special = mf.special_residues(consensusdict[prt], consensusdict[i])
		blosumcol = special["blosum80_score"]
		result.insert(dfsize - 1, i + "_blosum80", blosumcol)
	#count number of speciality occurrences for each position and add it to the result df
	columnnames = list(result.columns)
	#get only blosum score columns
	blosums = [x for x in columnnames if "blosum" in x]
	blosumdf = result[blosums]
	#we turn each row to a list, count number of blosum scores and add it to the df
	result = pd.concat([result, blosumdf.apply(lambda row: len([x for x in list(row) if x != "-"]), axis = 1)], axis = 1)
	result.rename(columns = {0 : "count"}, inplace = True)
	#print(result.tail(100))
	result.to_csv(currentcsvpath + prt + "_binary_comparisons.csv", index = False)
############################################################################################################################
print("My program took", time.time() - start_time, "seconds to run")
############################################################################################################################
