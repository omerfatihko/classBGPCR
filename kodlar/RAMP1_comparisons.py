from functions import my_functions as mf
import pandas as pd
import time
from os import listdir
from os.path import isfile, join
start_time = time.time()
pd.set_option("display.max_rows", None) #, "display.max_columns", None
"https://doi.org/10.1016/j.tips.2020.01.009 table 1 and https://doi.org/10.1016/j.apsb.2021.07.028 table 1 and 2 is used to assess RAMP1 interactions"

#combine CALCR, CALRL, GHRHR, GIPR, GCGR, GLP1R, GLP2R, PTH1R, PTH2R, PACR, VIPR1, VIPR2 and SCTR fasta files with orthologs
#It seems like almost all of them interact with RAMPs to some capacity, fuck
calcrfasta = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/calcr_orthologs_msa_pruned.fasta")
calrlfasta = mf.read_fasta("/cta/users/ofkonar/work/resources/class_B1/canonical/calrl_orthologs_msa_pruned.fasta")

RAMPfastas = calrlfasta + calcrfasta

#transform the list to a dataframe
RAMPdf = mf.fasta_to_dataframe(RAMPfastas)
#get the consensus for RAMP interacting GPCRs
consensus = mf.consensus(RAMPdf, 0.9)
#keep a deep copy so adjustments on the first one will not affect this one
rawconsensus = consensus.copy(deep = True)
#write the consensus down
#consensus.to_csv("/cta/users/ofkonar/work/results/csvs/RAMPinteracting_consensus.csv")

#make the comparisons
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
#add a prefix to the local consensus (RAMP interacting proteins consensus) so we can differentiate them
consensus = consensus.add_prefix("RAMP_")
#concatenate two consensus files, we will build upon this
result = pd.concat([globalconsensus, consensus], axis = 1)
#add human sequences to the data
#get the individual id and sequence names for CALCR and CALRL
calcrid = prttoid["CALCR"]
calrlid = prttoid["CALRL"]
calcrtitle = [x for x in humanseqs if calcrid in x]
calrltitle = [x for x in humanseqs if calrlid in x]
#get the sequences
calcrseq = tempdf[calcrtitle]
calrlseq = tempdf[calrltitle]
#they are one column dataframes so we squeeze them to series
calcrseq = calcrseq.squeeze()
calrlseq = calrlseq.squeeze()

#add calcr together with its sequence, wootten numbering and domain info
calcrcounter = 1
calcrWootten = []
calcrdomains = []
#add residue index to residue name, create Wootten column
#get the required protein column
calcrcolumn = Wootten_table["CALCR"]
for items in calcrseq.iteritems():
	if items[1] != "-":
		newres = items[1] + str(calcrcounter)
		calcrseq.set_value(items[0], newres)
		#get domain info of the residue
		a = domain_table.loc[domain_table["CALCR"] >=calcrcounter, "CALCR"]
		domain = a.index.values[0]
		calcrdomains.append(domain)
		calcrcounter += 1
		if newres in set(calcrcolumn):
			#get the index of Wootten
			windex = calcrcolumn[calcrcolumn == newres].index[0]
			#append to the Wootten column
			calcrWootten.append(str(Wootten_table.iloc[windex]['Wootten']))
		else:
			#if res lacks Wootten numbering
			calcrWootten.append('-')
	else:
		#for gaps
		calcrWootten.append('-')
		calcrdomains.append("-")
result.insert(result.shape[1], "CALCR_seq", calcrseq)
result.insert(result.shape[1], "CALCR_Wootten", calcrWootten, True)
result.insert(result.shape[1], "CALCR_Domain", calcrdomains, True)

#add calrl together with its sequence, wootten numbering and domain info
calrlcounter = 1
calrlWootten = []
calrldomains = []
#add residue index to residue name, create Wootten column
#get the required protein column
calrlcolumn = Wootten_table["CALRL"]
for items in calrlseq.iteritems():
	if items[1] != "-":
		newres = items[1] + str(calrlcounter)
		calrlseq.set_value(items[0], newres)
		#get domain info of the residue
		a = domain_table.loc[domain_table["CALRL"] >=calrlcounter, "CALRL"]
		domain = a.index.values[0]
		calrldomains.append(domain)
		calrlcounter += 1
		if newres in set(calrlcolumn):
			#get the index of Wootten
			windex = calrlcolumn[calrlcolumn == newres].index[0]
			#append to the Wootten column
			calrlWootten.append(str(Wootten_table.iloc[windex]['Wootten']))
		else:
			#if res lacks Wootten numbering
			calrlWootten.append('-')
	else:
		#for gaps
		calrlWootten.append('-')
		calrldomains.append("-")
result.insert(result.shape[1], "CALRL_seq", calrlseq)
result.insert(result.shape[1], "CALRL_Wootten", calrlWootten, True)
result.insert(result.shape[1], "CALRL_Domain", calrldomains, True)
#get the list of fasta files of proteins that we will compare with the RAMP interacting proteins
files = []
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/crfr1_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/crfr2_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/gcgr_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/ghrhr_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/gipr_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/glp1r_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/glp2r_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/pacr_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/pth1r_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/pth2r_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/sctr_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/vipr1_orthologs_msa_pruned.fasta")
files.append("/cta/users/ofkonar/work/resources/class_B1/canonical/vipr2_orthologs_msa_pruned.fasta")
#compare them one by one with the RAMP interacting protein consensus
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

result.to_csv("/cta/users/ofkonar/work/results/csvs/RAMPinteracting_binary_comparisons.csv", index = False)

print("My program took", time.time() - start_time, "seconds to run")