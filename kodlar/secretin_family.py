from functions import my_functions as mf
import pandas as pd
from os import listdir
from os.path import isfile, join
from collections import Counter
import re
unified_fasta_path = "/cta/users/ofkonar/work/resources/fasta/unified_fasta.fasta"
secretin_fasta_path = "/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/"
pd.set_option("display.max_rows", None, "display.max_columns", None)

#single operations
############################################################################################################################
#single orthologs aligned, and pruned

#ghrhr
#select orthologs
#ghrhrorthologs = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/ghrhr/ghrhr_orthologs.csv", 0)
#filter orthologs from the unified data
#ghrhrorthologsfastas = mf.filter_fasta(ghrhrorthologs, unified_fasta_path)
#write them down
#mf.write_fasta(secretin_fasta_path + "ghrhr_orthologs.fasta", ghrhrorthologsfastas)

#pacr
#select orthologs
#pacrorthologs = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/pacr/pacr_orthologs.csv", 0)
#filter orthologs from the unified data
#pacrorthologsfastas = mf.filter_fasta(pacrorthologs, unified_fasta_path)
#write them down
#mf.write_fasta(secretin_fasta_path + "pacr_orthologs.fasta", pacrorthologsfastas)

#sctr
#select orthologs
#sctrorthologs = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/sctr/sctr_orthologs.csv", 0)
#filter orthologs from the unified data
#sctrorthologsfastas = mf.filter_fasta(sctrorthologs, unified_fasta_path)
#write them down
#mf.write_fasta(secretin_fasta_path + "sctr_orthologs.fasta", sctrorthologsfastas)

#vipr1
#select orthologs
#vipr1orthologs = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/vipr1/vipr1_orthologs.csv", 0)
#filter orthologs from the unified data
#vipr1orthologsfastas = mf.filter_fasta(vipr1orthologs, unified_fasta_path)
#write them down
#mf.write_fasta(secretin_fasta_path + "vipr1_orthologs.fasta", vipr1orthologsfastas)

#vipr2
#select orthologs
#vipr2orthologs = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/vipr2/vipr2_orthologs.csv", 0)
#filter orthologs from the unified data
#vipr2orthologsfastas = mf.filter_fasta(vipr2orthologs, unified_fasta_path)
#write them down
#mf.write_fasta(secretin_fasta_path + "vipr2_orthologs.fasta", vipr2orthologsfastas)

############################################################################################################################
#prune the gaps and get consensus

#ghrhr
#read the aligned fasta file
#ghrhrorthologsfastas = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/ghrhr_orthologs_al.fasta")
#prune the gaps
#ghrhrorthologsfastaspruned = mf.fasta_pruner(ghrhrorthologsfastas)
#write down the pruned fasta for later use
#mf.write_fasta(secretin_fasta_path + "ghrhr_orthologs_al_pruned.fasta", ghrhrorthologsfastaspruned)
#turn the pruned fasta into dataframe since consensus function takes pandas dataframe
#ghrhrorthologsfastaspruneddf = mf.fasta_to_dataframe(ghrhrorthologsfastaspruned)
#get the consensus, there is no need for position adjustment since we already pruned gaps in the only human sequence
#ghrhrconsensus = mf.consensus(ghrhrorthologsfastaspruneddf, 0.9)
#write down the consensus sequence
#ghrhrconsensus.to_csv("/cta/users/ofkonar/work/results/secretin_subfamily/ghrhr/ghrhr_consensus.csv")

#pacr
#read the aligned fasta file
#pacrorthologsfastas = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pacr_orthologs_al.fasta")
#prune the gaps
#pacrorthologsfastaspruned = mf.fasta_pruner(pacrorthologsfastas)
#write down the pruned fasta for later use
#mf.write_fasta(secretin_fasta_path + "pacr_orthologs_al_pruned.fasta", pacrorthologsfastaspruned)
#turn the pruned fasta into dataframe since consensus function takes pandas dataframe
#pacrorthologsfastaspruneddf = mf.fasta_to_dataframe(pacrorthologsfastaspruned)
#get the consensus, there is no need for position adjustment since we already pruned gaps in the only human sequence
#pacrconsensus = mf.consensus(pacrorthologsfastaspruneddf, 0.9)
#write down the consensus sequence
#pacrconsensus.to_csv("/cta/users/ofkonar/work/results/secretin_subfamily/pacr/pacr_consensus.csv")

#sctr
#read the aligned fasta file
#sctrorthologsfastas = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/sctr_orthologs_al.fasta")
#prune the gaps
#sctrorthologsfastaspruned = mf.fasta_pruner(sctrorthologsfastas)
#write down the pruned fasta for later use
#mf.write_fasta(secretin_fasta_path + "sctr_orthologs_al_pruned.fasta", sctrorthologsfastaspruned)
#turn the pruned fasta into dataframe since consensus function takes pandas dataframe
#sctrorthologsfastaspruneddf = mf.fasta_to_dataframe(sctrorthologsfastaspruned)
#get the consensus, there is no need for position adjustment since we already pruned gaps in the only human sequence
#sctrconsensus = mf.consensus(sctrorthologsfastaspruneddf, 0.9)
#write down the consensus sequence
#sctrconsensus.to_csv("/cta/users/ofkonar/work/results/secretin_subfamily/sctr/sctr_consensus.csv")

#vipr1
#read the aligned fasta file
#vipr1orthologsfastas = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/vipr1_orthologs_al.fasta")
#prune the gaps
#vipr1orthologsfastaspruned = mf.fasta_pruner(vipr1orthologsfastas)
#write down the pruned fasta for later use
#mf.write_fasta(secretin_fasta_path + "vipr1_orthologs_al_pruned.fasta", vipr1orthologsfastaspruned)
#turn the pruned fasta into dataframe since consensus function takes pandas dataframe
#vipr1orthologsfastaspruneddf = mf.fasta_to_dataframe(vipr1orthologsfastaspruned)
#get the consensus, there is no need for position adjustment since we already pruned gaps in the only human sequence
#vipr1consensus = mf.consensus(vipr1orthologsfastaspruneddf, 0.9)
#write down the consensus sequence
#vipr1consensus.to_csv("/cta/users/ofkonar/work/results/secretin_subfamily/vipr1/vipr1_consensus.csv")

#vipr2
#read the aligned fasta file
#vipr2orthologsfastas = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/vipr2_orthologs_al.fasta")
#prune the gaps
#vipr2orthologsfastaspruned = mf.fasta_pruner(vipr2orthologsfastas)
#write down the pruned fasta for later use
#mf.write_fasta(secretin_fasta_path + "vipr2_orthologs_al_pruned.fasta", vipr2orthologsfastaspruned)
#turn the pruned fasta into dataframe since consensus function takes pandas dataframe
#vipr2orthologsfastaspruneddf = mf.fasta_to_dataframe(vipr2orthologsfastaspruned)
#get the consensus, there is no need for position adjustment since we already pruned gaps in the only human sequence
#vipr2consensus = mf.consensus(vipr2orthologsfastaspruneddf, 0.9)
#write down the consensus sequence
#vipr2consensus.to_csv("/cta/users/ofkonar/work/results/secretin_subfamily/vipr2/vipr2_consensus.csv")
############################################################################################################################
############################################################################################################################

#all orthologs together
############################################################################################################################
#combine all ortholog fastas for msa

#read fasta files
#ghrhrorthologs = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/ghrhr_orthologs.fasta")
#pacrorthologs = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pacr_orthologs.fasta")
#sctrorthologs = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/sctr_orthologs.fasta")
#vipr1orthologs = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/vipr1_orthologs.fasta")
#vipr2orthologs = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/vipr2_orthologs.fasta")
#allcombined = ghrhrorthologs + pacrorthologs + sctrorthologs + vipr1orthologs + vipr2orthologs
#mf.write_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_orthologs.fasta", allcombined)
############################################################################################################################
#prune aligned orthologs

#combinedorthologs = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_orthologs_al.fasta")

#get human sequences
#human_seqs = []
#for seq in combinedorthologs:
#	if "HUMAN" in seq:
#		human_seqs.append(seq)

#create a dictionary where header and body of fasta are split
#list_of_fastas = []
#for seq in combinedorthologs:
#	#get seperate lines
#	split_seq = seq.split("\n")
#	#get header (first line)
#	header = split_seq[0]
#	#combine the body
#	body = "".join(map(str,split_seq[1:]))
#	#assign header and body
#	list_of_fastas.append({"header": header, "body": body})

#create same dictionary for only human sequences
#list_of_humans = []
#for seq in human_seqs:
#	split_seq = seq.split("\n")
#	header = split_seq[0]
#	body = "".join(map(str,split_seq[1:]))
#	list_of_humans.append({"header": header, "body": body})

#get indices of gaps
#doublegapindices = []
#for i in range(len(list_of_humans[0]["body"])):
#	if (list_of_humans[0]["body"][i] == "-") & (list_of_humans[1]["body"][i] == "-") & (list_of_humans[2]["body"][i] == "-") & (list_of_humans[3]["body"][i] == "-") & (list_of_humans[4]["body"][i] == "-"):
#		doublegapindices.append(i)

#prune the body of fastas
#for i in range(len(list_of_fastas)):
	#get body
#	temp = list_of_fastas[i]["body"]
#	#remove double gap indices
#	pruned_body = mf.Pruner(temp, doublegapindices)
#	#reassign pruned body to body
#	list_of_fastas[i]["body"] = pruned_body

#return the fastas to list state
#fasta_like = []

#for element in list_of_fastas:
#	temp = element["header"] + "\n"
#	counter = (-(-len(element["body"])//60))-1
#	temp2 = element["body"]
#	for i in range(counter):
#		temp += (temp2[0:60] + "\n")
#		temp3 = temp2[60:]
#		temp2 = temp3
#	temp += (temp2 + "\n")
#	fasta_like.append(temp)
#mf.write_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_orthologs_al_pruned.fasta",fasta_like)
############################################################################################################################
#get aligned human sequences

#read aligned and pruned orthologs
#secretinorthologspruned = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_orthologs_al_pruned.fasta")
#get human sequences
#human_seqs = []
#for seq in secretinorthologspruned:
#	if "HUMAN" in seq:
#		human_seqs.append(seq)
#write them down
#mf.write_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_humans_aligned_pruned.fasta", human_seqs)
############################################################################################################################
#secretin family consensus

#read aligned and pruned orthologs
#secretinorthologspruned = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_orthologs_al_pruned.fasta")
#transform them into dataframe
#secretinorthologspruneddf = mf.fasta_to_dataframe(secretinorthologspruned)
#get consensus of secretin family
#secretinallconsensus = mf.consensus(secretinorthologspruneddf, 0.9)
#write them down
#secretinallconsensus.to_csv("/cta/users/ofkonar/work/results/secretin_subfamily/secretin_subfamily_orthologs_consensus.csv")

#get individual human protein orthologs from all aligned orthologs and calculate consensus with them, we will take the frequencies from them 
#ghrhr
#ghrhrorthologs = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/ghrhr/ghrhr_orthologs.csv", 0)
#ghrhrorthologsfastas = mf.filter_fasta(ghrhrorthologs, "/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_orthologs_al_pruned.fasta")
#ghrhrorthologsfastasdf = mf.fasta_to_dataframe(ghrhrorthologsfastas)
#ghrhrconsensus = mf.consensus(ghrhrorthologsfastasdf, 0.9)
#ghrhrfrequency = ghrhrconsensus["frequency"]
#pacr
#pacrorthologs = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/pacr/pacr_orthologs.csv", 0)
#pacrorthologsfastas = mf.filter_fasta(pacrorthologs, "/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_orthologs_al_pruned.fasta")
#pacrorthologsfastasdf = mf.fasta_to_dataframe(pacrorthologsfastas)
#pacrconsensus = mf.consensus(pacrorthologsfastasdf, 0.9)
#pacrfrequency = pacrconsensus["frequency"]
#sctr
#sctrorthologs = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/sctr/sctr_orthologs.csv", 0)
#sctrorthologsfastas = mf.filter_fasta(sctrorthologs, "/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_orthologs_al_pruned.fasta")
#sctrorthologsfastasdf = mf.fasta_to_dataframe(sctrorthologsfastas)
#sctrconsensus = mf.consensus(sctrorthologsfastasdf, 0.9)
#sctrfrequency = sctrconsensus["frequency"]
#vipr1
#vipr1orthologs = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/vipr1/vipr1_orthologs.csv", 0)
#vipr1orthologsfastas = mf.filter_fasta(vipr1orthologs, "/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_orthologs_al_pruned.fasta")
#vipr1orthologsfastasdf = mf.fasta_to_dataframe(vipr1orthologsfastas)
#vipr1consensus = mf.consensus(vipr1orthologsfastasdf, 0.9)
#vipr1frequency = vipr1consensus["frequency"]
#vipr2
#vipr2orthologs = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/vipr2/vipr2_orthologs.csv", 0)
#vipr2orthologsfastas = mf.filter_fasta(vipr2orthologs, "/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_orthologs_al_pruned.fasta")
#vipr2orthologsfastasdf = mf.fasta_to_dataframe(vipr2orthologsfastas)
#vipr2consensus = mf.consensus(vipr2orthologsfastasdf, 0.9)
#vipr2frequency = vipr2consensus["frequency"]

#take single human sequences and index them, we will their sequence and their index
#read aligned human sequences
#secretinhumansequences = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/secretin_subfamily_humans_aligned_pruned.fasta")
#get individual sequences
#human_ghrhr = []
#human_pacr = []
#human_sctr = []
#human_vipr1 = []
#human_vipr2 = []
#for seq in secretinhumansequences:
#	if "GHRHR" in seq:
#		human_ghrhr.append(seq)
#	elif "PACR" in seq:
#		human_pacr.append(seq)
#	elif "SCTR" in seq:
#		human_sctr.append(seq)
#	elif "VIPR1" in seq:
#		human_vipr1.append(seq)
#	elif "VIPR2" in seq:
#		human_vipr2.append(seq)

#turn them into dataframes
#human_ghrhrdf = mf.fasta_to_dataframe(human_ghrhr)
#human_pacrdf = mf.fasta_to_dataframe(human_pacr)
#human_sctrdf = mf.fasta_to_dataframe(human_sctr)
#human_vipr1df = mf.fasta_to_dataframe(human_vipr1)
#human_vipr2df = mf.fasta_to_dataframe(human_vipr2)

#index the sequences
#ghrhr
#index = 1
#ghrhrindex = []
#for i in range(len(human_ghrhrdf)):
#	if human_ghrhrdf.iloc[i,0] == "-":
#		ghrhrindex.append("-")
#	else:
#		ghrhrindex.append(index)
#		index+=1
#add index and frequency columns and rename the first column to the protein name
#human_ghrhrdf["GHRHR_index"] = ghrhrindex
#human_ghrhrdf["GHRHR_frequency"]= ghrhrfrequency
#human_ghrhrdf.rename(columns = {human_ghrhrdf.columns[0]: "GHRHR"}, inplace = True)
#pacr
#index = 1
#pacrindex = []
#for i in range(len(human_pacrdf)):
#	if human_pacrdf.iloc[i,0] == "-":
#		pacrindex.append("-")
#	else:
#		pacrindex.append(index)
#		index+=1
#add index and frequency columns and rename the first column to the protein name
#human_pacrdf["PACR_index"] = pacrindex
#human_pacrdf["PACR_frequency"]= pacrfrequency
#human_pacrdf.rename(columns = {human_pacrdf.columns[0]: "PACR"}, inplace = True)
#sctr
#index = 1
#sctrindex = []
#for i in range(len(human_sctrdf)):
#	if human_sctrdf.iloc[i,0] == "-":
#		sctrindex.append("-")
#	else:
#		sctrindex.append(index)
#		index+=1
#add index and frequency columns and rename the first column to the protein name
#human_sctrdf["SCTR_index"] = sctrindex
#human_sctrdf["SCTR_frequency"]= sctrfrequency
#human_sctrdf.rename(columns = {human_sctrdf.columns[0]: "SCTR"}, inplace = True)
#vipr1
#index = 1
#vipr1index = []
#for i in range(len(human_vipr1df)):
#	if human_vipr1df.iloc[i,0] == "-":
#		vipr1index.append("-")
#	else:
#		vipr1index.append(index)
#		index+=1
#add index and frequency columns and rename the first column to the protein name
#human_vipr1df["VIPR1_index"] = vipr1index
#human_vipr1df["VIPR1_frequency"]= vipr1frequency
#human_vipr1df.rename(columns = {human_vipr1df.columns[0]: "VIPR1"}, inplace = True)
#vipr2
#index = 1
#vipr2index = []
#for i in range(len(human_vipr2df)):
#	if human_vipr2df.iloc[i,0] == "-":
#		vipr2index.append("-")
#	else:
#		vipr2index.append(index)
#		index+=1
#add index and frequency columns and rename the first column to the protein name
#human_vipr2df["VIPR2_index"] = vipr2index
#human_vipr2df["VIPR2_frequency"]= vipr2frequency
#human_vipr2df.rename(columns = {human_vipr2df.columns[0]: "VIPR2"}, inplace = True)

#combine individual consensus sequences with secretin subfamily consensus sequence
#combined = pd.concat([human_ghrhrdf, human_pacrdf, human_sctrdf, human_vipr1df, human_vipr2df, secretinallconsensus], axis = 1)
#combined.to_csv("/cta/users/ofkonar/work/results/secretin_subfamily/secretin_subfamily_orthologs_consensus_withhuman.csv")
############################################################################################################################
############################################################################################################################

#within secretin family pair comparisons
############################################################################################################################
#combine all pairs
#read orthologs
#ghrhr = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/ghrhr_orthologs.fasta")
#pacr = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pacr_orthologs.fasta")
#sctr = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/sctr_orthologs.fasta")
#vipr1 = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/vipr1_orthologs.fasta")
#vipr2 = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/vipr2_orthologs.fasta")
#form a list of lists of orthologs
#all_list = [ghrhr, pacr, sctr, vipr1, vipr2]
#specify write path
#path = "/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/"
#form a list of protein names
#prt_names = ["ghrhr", "pacr", "sctr", "vipr1", "vipr2"]
#combine every list of orthologs with every other list of orthologs
#for i in range(len(all_list)):
#	for j in range(i+1, len(all_list)):
#		#combine ortholog lists
#		temp = all_list[i] + all_list[j]
#		#define exact write path
#		write_path = path + prt_names[i] + "_" + prt_names[j] + "_orthologs_combined.fasta"
#		#write them down
#		mf.write_fasta( write_path, temp)
############################################################################################################################
#"/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/ghrhr_pacr_orthologs_combined.fasta"
#"/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/ghrhr_sctr_orthologs_combined.fasta"
#"/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/ghrhr_vipr1_orthologs_combined.fasta"
#"/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/ghrhr_vipr2_orthologs_combined.fasta"
#"/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/pacr_sctr_orthologs_combined.fasta"
#"/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/pacr_vipr1_orthologs_combined.fasta"
#"/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/pacr_vipr2_orthologs_combined.fasta"
#"/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/sctr_vipr1_orthologs_combined.fasta"
#"/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/sctr_vipr2_orthologs_combined.fasta"
#"/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/vipr1_vipr2_orthologs_combined.fasta"
############################################################################################################################
#split and compare the pairs
#pair_path = "/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pairs/"
results_path = "/cta/users/ofkonar/work/results/secretin_subfamily/"

#get all aligned ortholog pairs
#fastafiles = [f for f in listdir(pair_path) if isfile(join(pair_path, f))]
#alignedfiles = [i for i in fastafiles if "_aligned" in i]

#for file in alignedfiles:
#	#read aligned fasta
#	filepath = pair_path + file
#	fasta = mf.read_fasta(filepath)
#	#prune aligned fasta
#	prunedfasta = mf.fasta_pruner(fasta)
#	#get protein names
#	parts = file.split("_")
#	prt1 = parts[0]
#	prt2 = parts[1]
#	#get ortholog list paths
#	orthologlist1 = results_path + prt1 + "/" + prt1 + "_orthologs.csv"
#	orthologlist2 = results_path + prt2 + "/" + prt2 + "_orthologs.csv"
#	#split the fasta into two parts
#	splitdictionary = mf.fasta_splitter(prunedfasta, orthologlist1, orthologlist2)
#	#get the parts
#	fastalist1 = splitdictionary["fasta_list1"]
#	fastalist2 = splitdictionary["fasta_list2"]
#	#turn the lists into dataframes
#	fastadf1 = mf.fasta_to_dataframe(fastalist1)
#	fastadf2 = mf.fasta_to_dataframe(fastalist2)
#	#get the consensus
#	consensus1 = mf.consensus(fastadf1, 0.9)
#	consensus2 = mf.consensus(fastadf2, 0.9)
#	#compare them
#	special1 = mf.special_residues(consensus1, consensus2)
#	special2 = mf.special_residues(consensus2, consensus1)
#	#finalize the results
#	final1 = mf.finalize(fastadf1, special1, prt1.upper())
#	final2 = mf.finalize(fastadf2, special2, prt2.upper())
#	#write the results down
#	writepath1 = results_path + prt1 + "/" + prt1 + "_vs_" + prt2 + "_special_residues.csv"
#	writepath2 = results_path + prt2 + "/" + prt2 + "_vs_" + prt1 + "_special_residues.csv"
#	final1.to_csv(writepath1)
#	final2.to_csv(writepath2)
############################################################################################################################
#compile the results and get quadruple special residues and others (triple, double, single)
#ghrhr
#ghrhr_path = "/cta/users/ofkonar/work/results/secretin_subfamily/ghrhr/"
#get result files, collect them in a dictionary
#onlyfiles = [f for f in listdir(ghrhr_path) if isfile(join(ghrhr_path, f))]
#ghrhr_result_files = [i for i in onlyfiles if "special_" in i]
#ghrhr_results = {}
#for file in ghrhr_result_files:
#	filepath = ghrhr_path + file
#	result = pd.read_csv(filepath, header = 0, index_col = 0)
#	parts = file.split("_")
#	companion = parts[2]
#	ghrhr_results[companion] = result

#count the occurance of residues
#get the names as lists and concatanate them
#residue_names = []
#for key in ghrhr_results:
#	value = ghrhr_results[key]
#	temp = value["Name"].to_list()
#	residue_names.extend(temp)
#Use counter to count occurences of residues
#counts = Counter(residue_names)
#get quadruples
#quadruples = []
#for key in counts:
#	if counts[key] == 4:
#		quadruples.append(key)
#collect the results
#quadruple_results = pd.DataFrame(columns = ["residue", "frequency", "pacr_score", "sctr_score", "vipr1_score", "vipr2_score","Wootten", "name", "domain"])
#temp = re.compile("([a-zA-Z]+)([0-9]+)")
#for element in quadruples:
#	parts = temp.match(element).groups()
#	residue = parts[0]
#	name = element
#	for key in ghrhr_results:
#		value = ghrhr_results[key]
#		if key == "pacr":
#			frequency = value.loc[value["Name"] == element, "frequency"].values[0]
#			pacr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#			Wootten = value.loc[value["Name"] == element, "Wootten"].values[0]
#			domain = value.loc[value["Name"] == element, "domains"].values[0]
#		if key == "sctr":
#			sctr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#		if key == "vipr1":
#			vipr1_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#		if key == "vipr2":
#			vipr2_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#	quadruple_results = quadruple_results.append({"residue" : residue, "frequency" : frequency, "pacr_score" : pacr_score, "sctr_score" : sctr_score, "vipr1_score" : vipr1_score, "vipr2_score" : vipr2_score, "Wootten" : Wootten, "name" : name, "domain" : domain}, ignore_index = True)
#quadruple_results.to_csv(ghrhr_path + "ghrhr_quadruple_special.csv")

#pacr
#pacr_path = "/cta/users/ofkonar/work/results/secretin_subfamily/pacr/"
#get result files, collect them in a dictionary
#onlyfiles = [f for f in listdir(pacr_path) if isfile(join(pacr_path, f))]
#pacr_result_files = [i for i in onlyfiles if "special_" in i]
#pacr_results = {}
#for file in pacr_result_files:
#	filepath = pacr_path + file
#	result = pd.read_csv(filepath, header = 0, index_col = 0)
#	parts = file.split("_")
#	companion = parts[2]
#	pacr_results[companion] = result

#count the occurance of residues
#get the names as lists and concatanate them
#residue_names = []
#for key in pacr_results:
#	value = pacr_results[key]
#	temp = value["Name"].to_list()
#	residue_names.extend(temp)
#Use counter to count occurences of residues
#counts = Counter(residue_names)
#print(counts)
#get quadruples
#quadruples = []
#for key in counts:
#	if counts[key] == 4:
#		quadruples.append(key)
#collect the results
#quadruple_results = pd.DataFrame(columns = ["residue", "frequency", "ghrhr_score", "sctr_score", "vipr1_score", "vipr2_score","Wootten", "name", "domain"])
#temp = re.compile("([a-zA-Z]+)([0-9]+)")
#for element in quadruples:
#	parts = temp.match(element).groups()
#	residue = parts[0]
#	name = element
#	for key in pacr_results:
#		value = pacr_results[key]
#		if key == "ghrhr":
#			frequency = value.loc[value["Name"] == element, "frequency"].values[0]
#			ghrhr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#			Wootten = value.loc[value["Name"] == element, "Wootten"].values[0]
#			domain = value.loc[value["Name"] == element, "domains"].values[0]
#		if key == "sctr":
#			sctr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#		if key == "vipr1":
#			vipr1_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#		if key == "vipr2":
#			vipr2_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#	quadruple_results = quadruple_results.append({"residue" : residue, "frequency" : frequency, "ghrhr_score" : ghrhr_score, "sctr_score" : sctr_score, "vipr1_score" : vipr1_score, "vipr2_score" : vipr2_score, "Wootten" : Wootten, "name" : name, "domain" : domain}, ignore_index = True)
#quadruple_results.to_csv(pacr_path + "pacr_quadruple_special.csv")

#sctr
#sctr_path = "/cta/users/ofkonar/work/results/secretin_subfamily/sctr/"
#get result files, collect them in a dictionary
#onlyfiles = [f for f in listdir(sctr_path) if isfile(join(sctr_path, f))]
#sctr_result_files = [i for i in onlyfiles if "special_" in i]
#sctr_results = {}
#for file in sctr_result_files:
#	filepath = sctr_path + file
#	result = pd.read_csv(filepath, header = 0, index_col = 0)
#	parts = file.split("_")
#	companion = parts[2]
#	sctr_results[companion] = result

#count the occurance of residues
#get the names as lists and concatanate them
#residue_names = []
#for key in sctr_results:
#	value = sctr_results[key]
#	temp = value["Name"].to_list()
#	residue_names.extend(temp)
#Use counter to count occurences of residues
#counts = Counter(residue_names)
#print(counts)
#get quadruples
#quadruples = []
#for key in counts:
#	if counts[key] == 4:
#		quadruples.append(key)
#collect the results
#quadruple_results = pd.DataFrame(columns = ["residue", "frequency", "ghrhr_score", "pacr_score", "vipr1_score", "vipr2_score","Wootten", "name", "domain"])
#temp = re.compile("([a-zA-Z]+)([0-9]+)")
#for element in quadruples:
#	parts = temp.match(element).groups()
#	residue = parts[0]
#	name = element
#	for key in sctr_results:
#		value = sctr_results[key]
#		if key == "ghrhr":
#			frequency = value.loc[value["Name"] == element, "frequency"].values[0]
#			ghrhr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#			Wootten = value.loc[value["Name"] == element, "Wootten"].values[0]
#			domain = value.loc[value["Name"] == element, "domains"].values[0]
#		if key == "pacr":
#			pacr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#		if key == "vipr1":
#			vipr1_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#		if key == "vipr2":
#			vipr2_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#	quadruple_results = quadruple_results.append({"residue" : residue, "frequency" : frequency, "ghrhr_score" : ghrhr_score, "pacr_score" : pacr_score, "vipr1_score" : vipr1_score, "vipr2_score" : vipr2_score, "Wootten" : Wootten, "name" : name, "domain" : domain}, ignore_index = True)
#quadruple_results.to_csv(sctr_path + "sctr_quadruple_special.csv")

#vipr1
#vipr1_path = "/cta/users/ofkonar/work/results/secretin_subfamily/vipr1/"
#get result files, collect them in a dictionary
#onlyfiles = [f for f in listdir(vipr1_path) if isfile(join(vipr1_path, f))]
#vipr1_result_files = [i for i in onlyfiles if "special_" in i]
#vipr1_results = {}
#for file in vipr1_result_files:
#	filepath = vipr1_path + file
#	result = pd.read_csv(filepath, header = 0, index_col = 0)
#	parts = file.split("_")
#	companion = parts[2]
#	vipr1_results[companion] = result

#count the occurance of residues
#get the names as lists and concatanate them
#residue_names = []
#for key in vipr1_results:
#	value = vipr1_results[key]
#	temp = value["Name"].to_list()
#	residue_names.extend(temp)
#Use counter to count occurences of residues
#counts = Counter(residue_names)
#print(counts)
#get quadruples
#quadruples = []
#for key in counts:
#	if counts[key] == 4:
#		quadruples.append(key)
#collect the results
#quadruple_results = pd.DataFrame(columns = ["residue", "frequency", "ghrhr_score", "pacr_score", "sctr_score", "vipr2_score","Wootten", "name", "domain"])
#temp = re.compile("([a-zA-Z]+)([0-9]+)")
#for element in quadruples:
#	parts = temp.match(element).groups()
#	residue = parts[0]
#	name = element
#	for key in vipr1_results:
#		value = vipr1_results[key]
#		if key == "ghrhr":
#			frequency = value.loc[value["Name"] == element, "frequency"].values[0]
#			ghrhr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#			Wootten = value.loc[value["Name"] == element, "Wootten"].values[0]
#			domain = value.loc[value["Name"] == element, "domains"].values[0]
#		if key == "pacr":
#			pacr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#		if key == "sctr":
#			sctr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#		if key == "vipr2":
#			vipr2_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#	quadruple_results = quadruple_results.append({"residue" : residue, "frequency" : frequency, "ghrhr_score" : ghrhr_score, "pacr_score" : pacr_score, "sctr_score" : sctr_score, "vipr2_score" : vipr2_score, "Wootten" : Wootten, "name" : name, "domain" : domain}, ignore_index = True)
#quadruple_results.to_csv(vipr1_path + "vipr1_quadruple_special.csv")

#vipr2
#vipr2_path = "/cta/users/ofkonar/work/results/secretin_subfamily/vipr2/"
#get result files, collect them in a dictionary
#onlyfiles = [f for f in listdir(vipr2_path) if isfile(join(vipr2_path, f))]
#vipr2_result_files = [i for i in onlyfiles if "special_" in i]
#vipr2_results = {}
#for file in vipr2_result_files:
#	filepath = vipr2_path + file
#	result = pd.read_csv(filepath, header = 0, index_col = 0)
#	parts = file.split("_")
#	companion = parts[2]
#	vipr2_results[companion] = result

#count the occurance of residues
#get the names as lists and concatanate them
#residue_names = []
#for key in vipr2_results:
#	value = vipr2_results[key]
#	temp = value["Name"].to_list()
#	residue_names.extend(temp)
#Use counter to count occurences of residues
#counts = Counter(residue_names)
#print(counts)
#get quadruples
#quadruples = []
#for key in counts:
#	if counts[key] == 4:
#		quadruples.append(key)
#collect the results
#quadruple_results = pd.DataFrame(columns = ["residue", "frequency", "ghrhr_score", "pacr_score", "sctr_score", "vipr1_score","Wootten", "name", "domain"])
#temp = re.compile("([a-zA-Z]+)([0-9]+)")
#for element in quadruples:
#	parts = temp.match(element).groups()
#	residue = parts[0]
#	name = element
#	for key in vipr2_results:
#		value = vipr2_results[key]
#		if key == "ghrhr":
#			frequency = value.loc[value["Name"] == element, "frequency"].values[0]
#			ghrhr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#			Wootten = value.loc[value["Name"] == element, "Wootten"].values[0]
#			domain = value.loc[value["Name"] == element, "domains"].values[0]
#		if key == "pacr":
#			pacr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#		if key == "sctr":
#			sctr_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#		if key == "vipr1":
#			vipr1_score = value.loc[value["Name"] == element, "blosum80_score"].values[0]
#	quadruple_results = quadruple_results.append({"residue" : residue, "frequency" : frequency, "ghrhr_score" : ghrhr_score, "pacr_score" : pacr_score, "sctr_score" : sctr_score, "vipr1_score" : vipr1_score, "Wootten" : Wootten, "name" : name, "domain" : domain}, ignore_index = True)
#quadruple_results.to_csv(vipr2_path + "vipr2_quadruple_special.csv")
############################################################################################################################
#special residues to look for
#ghrhr
#path to get quadruple specials
ghrhr_path = results_path + "ghrhr/"
#get files in the path
ghrhrfiles =  [f for f in listdir(ghrhr_path) if isfile(join(ghrhr_path, f))]
#get file with quadruple special residues
ghrhrquadfile = [i for i in ghrhrfiles if "quadruple" in i]
#read the file
ghrhrquadresults = pd.read_csv(ghrhr_path + ghrhrquadfile[0], header = 0, index_col = 0)
#get residue names
ghrhrnames = ghrhrquadresults["name"].to_list()
#we need residue indexes so we seperate letter code and index using re
ghrhrlookfor = []
temp = re.compile("([a-zA-Z]+)([0-9]+)")
for element in ghrhrnames:
	parts = temp.match(element).groups()
	index = int(parts[1])
	ghrhrlookfor.append(index)
#get the results with lookout
ghrhrresults= mf.lookout("secretin", "GHRHR", ghrhrlookfor)
ghrhrwritepath = ghrhr_path + "ghrhr_lookout.csv"
ghrhrresults.to_csv(ghrhrwritepath)

#pacr
#path to get quadruple specials
pacr_path = results_path + "pacr/"
#get files in the path
pacrfiles =  [f for f in listdir(pacr_path) if isfile(join(pacr_path, f))]
#get file with quadruple special residues
pacrquadfile = [i for i in pacrfiles if "quadruple" in i]
#read the file
pacrquadresults = pd.read_csv(pacr_path + pacrquadfile[0], header = 0, index_col = 0)
#get residue names
pacrnames = pacrquadresults["name"].to_list()
#we need residue indexes so we seperate letter code and index using re
pacrlookfor = []
temp = re.compile("([a-zA-Z]+)([0-9]+)")
for element in pacrnames:
	parts = temp.match(element).groups()
	index = int(parts[1])
	pacrlookfor.append(index)
#get the results with lookout
pacrresults= mf.lookout("secretin", "PACR", pacrlookfor)
pacrwritepath = pacr_path + "pacr_lookout.csv"
pacrresults.to_csv(pacrwritepath)

#sctr
#path to get quadruple specials
sctr_path = results_path + "sctr/"
#get files in the path
sctrfiles =  [f for f in listdir(sctr_path) if isfile(join(sctr_path, f))]
#get file with quadruple special residues
sctrquadfile = [i for i in sctrfiles if "quadruple" in i]
#read the file
sctrquadresults = pd.read_csv(sctr_path + sctrquadfile[0], header = 0, index_col = 0)
#get residue names
sctrnames = sctrquadresults["name"].to_list()
#we need residue indexes so we seperate letter code and index using re
sctrlookfor = []
temp = re.compile("([a-zA-Z]+)([0-9]+)")
for element in sctrnames:
	parts = temp.match(element).groups()
	index = int(parts[1])
	sctrlookfor.append(index)
#get the results with lookout
sctrresults= mf.lookout("secretin", "SCTR", sctrlookfor)
sctrwritepath = sctr_path + "sctr_lookout.csv"
sctrresults.to_csv(sctrwritepath)

#vipr1
#path to get quadruple specials
vipr1_path = results_path + "vipr1/"
#get files in the path
vipr1files =  [f for f in listdir(vipr1_path) if isfile(join(vipr1_path, f))]
#get file with quadruple special residues
vipr1quadfile = [i for i in vipr1files if "quadruple" in i]
#read the file
vipr1quadresults = pd.read_csv(vipr1_path + vipr1quadfile[0], header = 0, index_col = 0)
#get residue names
vipr1names = vipr1quadresults["name"].to_list()
#we need residue indexes so we seperate letter code and index using re
vipr1lookfor = []
temp = re.compile("([a-zA-Z]+)([0-9]+)")
for element in vipr1names:
	parts = temp.match(element).groups()
	index = int(parts[1])
	vipr1lookfor.append(index)
#get the results with lookout
vipr1results= mf.lookout("secretin", "VIPR1", vipr1lookfor)
vipr1writepath = vipr1_path + "vipr1_lookout.csv"
vipr1results.to_csv(vipr1writepath)

#vipr2
#path to get quadruple specials
vipr2_path = results_path + "vipr2/"
#get files in the path
vipr2files =  [f for f in listdir(vipr2_path) if isfile(join(vipr2_path, f))]
#get file with quadruple special residues
vipr2quadfile = [i for i in vipr2files if "quadruple" in i]
#read the file
vipr2quadresults = pd.read_csv(vipr2_path + vipr2quadfile[0], header = 0, index_col = 0)
#get residue names
vipr2names = vipr2quadresults["name"].to_list()
#we need residue indexes so we seperate letter code and index using re
vipr2lookfor = []
temp = re.compile("([a-zA-Z]+)([0-9]+)")
for element in vipr2names:
	parts = temp.match(element).groups()
	index = int(parts[1])
	vipr2lookfor.append(index)
#get the results with lookout
vipr2results= mf.lookout("secretin", "VIPR2", vipr2lookfor)
vipr2writepath = vipr2_path + "vipr2_lookout.csv"
vipr2results.to_csv(vipr2writepath)