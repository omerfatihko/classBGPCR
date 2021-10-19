from functions import my_functions as mf
import pandas as pd
import time
from os import listdir
from os.path import isfile, join
start_time = time.time()
#according to STRING database, GCGR (GLR), GIPR, GLP1R, GLP2R, PACR, SCTR, VIPR1, and VIPR2 interact with GCG (glucagon)
#Although GHRHR is part of GCGR like family it does not.
gcglist = ["GCGR", "GIPR", "GLP1R", "GLP2R", "PACR", "SCTR", "VIPR1", "VIPR2"]
#combine afformentioned GPCRs ortholog files
gcgrfasta = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/gcgr_orthologs_msa_pruned.fasta")
giprfasta = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/gipr_orthologs_msa_pruned.fasta")
glp1rfasta = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/glp1r_orthologs_msa_pruned.fasta")
glp2rfasta = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/glp2r_orthologs_msa_pruned.fasta")
pacrfasta = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/pacr_orthologs_msa_pruned.fasta")
sctrfasta = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/sctr_orthologs_msa_pruned.fasta")
vipr1fasta = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/vipr1_orthologs_msa_pruned.fasta")
vipr2fasta = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/vipr2_orthologs_msa_pruned.fasta")
gcginteractingorthologs = gcgrfasta + giprfasta + glp1rfasta + glp2rfasta + pacrfasta + sctrfasta + vipr1fasta + vipr2fasta

#transform the list to a dataframe
gcginteractingdf = mf.fasta_to_dataframe(gcginteractingorthologs)
#get the consensus for GCG interacting GPCRs
consensus = mf.consensus(gcginteractingdf, 0.9)
#keep a deep copy so adjustments on the first one will not affect this one
rawconsensus = consensus.copy(deep = True)
#write the consensus down
#consensus.to_csv("/cta/users/ofkonar/work/results/csvs/GCGinteracting_consensus.csv")

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

#read the Wootten file for later use
Wootten_table = pd.read_csv("/cta/users/ofkonar/work/resources/residue_table.csv")
#read the domain file for later use
domain_table = pd.read_csv("/cta/users/ofkonar/work/resources/class_B1_domains.csv", index_col = 0)
#read the global consensus 
globalconsensus = pd.read_csv("/cta/users/ofkonar/work/results/csvs/class_B1_global_consensus.csv", index_col = 0)
#add a prefix to the local consensus (GCG interacting proteins consensus) so we can differentiate them
consensus = consensus.add_prefix("GCGi_")
#concatenate two consensus files, we will build upon this
result = pd.concat([globalconsensus, consensus], axis = 1)

#add human sequences to the data
#get the individual id and sequence names for GCG interacting GPCRs: GCGR (GLR), GIPR, GLP1R, GLP2R, PACR, SCTR, VIPR1, and VIPR2
for i in gcglist:
	interactorid = prttoid[i]
	interactortitle = [x for x in humanseqs if interactorid in x]
	#get the sequences
	interactorseq = tempdf[interactortitle]
	#they are one column dataframes so we squeeze them to series
	interactorseq = interactorseq.squeeze()

	#add interactor together with its sequence, wootten numbering and domain info
	counter = 1
	Wootten = []
	domains = []
	#add residue index to residue name, create Wootten column
	#get the required protein column
	interactorcolumn = Wootten_table[i]
	for items in interactorseq.iteritems():
		if items[1] != "-":
			newres = items[1] + str(counter)
			interactorseq.set_value(items[0], newres)
			#get domain info of the residue
			a = domain_table.loc[domain_table[i] >=counter, i]
			domain = a.index.values[0]
			domains.append(domain)
			counter += 1
			if newres in set(interactorcolumn):
				#get the index of Wootten
				windex = interactorcolumn[interactorcolumn == newres].index[0]
				#append to the Wootten column
				Wootten.append(str(Wootten_table.iloc[windex]['Wootten']))
			else:
				#if res lacks Wootten numbering
				Wootten.append('-')
		else:
			#for gaps
			Wootten.append('-')
			domains.append("-")
	result.insert(result.shape[1], i + "_seq", interactorseq)
	result.insert(result.shape[1], i + "_Wootten", Wootten, True)
	result.insert(result.shape[1], i + "_Domain", domains, True)

#get the list of fasta files of proteins that we will compare with the RAMP interacting proteins: CALCR, CALRL, CRFR1, CRFR2, GHRHR, PTH1R, PTH2R
files = ["/cta/users/ofkonar/work/resources/class_B1/canonical/calcr_orthologs_msa_pruned.fasta",
"/cta/users/ofkonar/work/resources/class_B1/canonical/calrl_orthologs_msa_pruned.fasta",
"/cta/users/ofkonar/work/resources/class_B1/canonical/crfr1_orthologs_msa_pruned.fasta",
"/cta/users/ofkonar/work/resources/class_B1/canonical/crfr2_orthologs_msa_pruned.fasta",
"/cta/users/ofkonar/work/resources/class_B1/canonical/ghrhr_orthologs_msa_pruned.fasta",
"/cta/users/ofkonar/work/resources/class_B1/canonical/pth1r_orthologs_msa_pruned.fasta",
"/cta/users/ofkonar/work/resources/class_B1/canonical/pth2r_orthologs_msa_pruned.fasta"]

#compare them one by one with the GCG interacting protein consensus
for i in files:
	#get compared proteins name
	prtname = i.split("/")[8].split("_")[0]
	#read the fasta file
	tempfasta = mf.read_fasta(i)
	#get the consensus of the protein
	fastadf = mf.fasta_to_dataframe(tempfasta)
	localconsensus = mf.consensus(fastadf, 0.9)
	#compare them
	special = mf.special_residues(rawconsensus, localconsensus)
	#get the blosum scores
	blosumcol = special["blosum80_score"]
	#insert it to the table
	result.insert(result.shape[1], prtname + "_blosum80", blosumcol)

#count number of speciality occurrences for each position and add it to the result df
columnnames = list(result.columns)
#get only blosum score columns
blosums = [x for x in columnnames if "blosum" in x]
blosumdf = result[blosums]
#we turn each row to a list, count number of blosum scores and add it to the df
result = pd.concat([result, blosumdf.apply(lambda row: len([x for x in list(row) if x != "-"]), axis = 1)], axis = 1)
result.rename(columns = {0 : "count"}, inplace = True)
print(result.head())
result.to_csv("/cta/users/ofkonar/work/results/csvs/GCGinteracting_binary_comparisons.csv", index = False)

print("My program took", time.time() - start_time, "seconds to run")