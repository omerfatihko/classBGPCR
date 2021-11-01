import pandas as pd
from functions import my_functions as mf

unified_fasta_path = "/cta/users/ofkonar/work/resources/fasta/unified_fasta.fasta"

#all_orthologs = mf.select_csv("/cta/users/ofkonar/work/results/all_orthologs.csv", 0)
#all_orthologs_fastas = mf.filter_fasta(all_orthologs, unified_fasta_path)
#mf.write_fasta("/cta/users/ofkonar/work/resources/fasta/all_orthologs.fasta", all_orthologs_fastas)
"""
all_orthologs = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/msa_all_orthologs.fasta")

#get human sequences
human_seqs = []
for seq in all_orthologs:
	if "HUMAN" in seq:
		human_seqs.append(seq)

#create a dictionary where header and body of fasta are split
list_of_fastas = []
for seq in all_orthologs:
	#get seperate lines
	split_seq = seq.split("\n")
	#get header (first line)
	header = split_seq[0]
	#combine the body
	body = "".join(map(str,split_seq[1:]))
	#assign header and body
	list_of_fastas.append({"header": header, "body": body})

#create same dictionary for only human sequences
list_of_humans = []
for seq in human_seqs:
	split_seq = seq.split("\n")
	header = split_seq[0]
	body = "".join(map(str,split_seq[1:]))
	list_of_humans.append({"header": header, "body": body})

#get indices of gaps
doublegapindices = []
for i in range(len(list_of_humans[0]["body"])):
	if (list_of_humans[0]["body"][i] == "-") & (list_of_humans[1]["body"][i] == "-") & (list_of_humans[2]["body"][i] == "-") & (list_of_humans[3]["body"][i] == "-"):
		doublegapindices.append(i)

#prune the body of fastas
for i in range(len(list_of_fastas)):
	#get body
	temp = list_of_fastas[i]["body"]
	#remove double gap indices
	pruned_body = mf.Pruner(temp, doublegapindices)
	#reassign pruned body to body
	list_of_fastas[i]["body"] = pruned_body

#return the fastas to list state
fasta_like = []

for element in list_of_fastas:
	temp = element["header"] + "\n"
	counter = (-(-len(element["body"])//60))-1
	temp2 = element["body"]
	for i in range(counter):
		temp += (temp2[0:60] + "\n")
		temp3 = temp2[60:]
		temp2 = temp3
	temp += (temp2 + "\n")
	fasta_like.append(temp)
"""
#mf.write_fasta("/cta/users/ofkonar/work/resources/fasta/all_orthologs_pruned.fasta",fasta_like)

#all_orthologs_pruned = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/all_orthologs_pruned.fasta")
"""
#get human sequences
human_seqs = []
for seq in all_orthologs_pruned:
	if "HUMAN" in seq:
		human_seqs.append(seq)

mf.write_fasta("/cta/users/ofkonar/work/resources/fasta/all_human_aligned.fasta", human_seqs)
"""
#aopdf = mf.fasta_to_dataframe(all_orthologs_pruned)
#aopdfc = mf.consensus(aopdf, 0.9)
#aopdfc.to_csv("/cta/users/ofkonar/work/results/all_orthologs_pruned_consensus.csv")

#all_human_aligned = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/all_human_aligned.fasta")
#ahadf = mf.fasta_to_dataframe(all_human_aligned)

#combined = pd.concat([ahadf, aopdfc], axis = 1)
#combined.to_csv("/cta/users/ofkonar/work/results/all_orthologs_pruned_consensus_withhuman.csv")

"""
gcgr_orthologs = mf.select_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_orthologs.csv", 0)
gcgr_orthologs_fastas = mf.filter_fasta(gcgr_orthologs, "/cta/users/ofkonar/work/resources/fasta/all_orthologs_pruned.fasta")
gcgr_orthologs_fastas_df = mf.fasta_to_dataframe(gcgr_orthologs_fastas)
gcgr_consensus = mf.consensus(gcgr_orthologs_fastas_df, 0.9)
gcgr_consensus_posadjusted = mf.consensus_posadjusted(gcgr_orthologs_fastas_df, gcgr_consensus)
gcgr_consensus.to_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_consensus.csv")
gcgr_consensus_posadjusted.to_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_consensus_posadjusted.csv")

gipr_orthologs = mf.select_csv("/cta/users/ofkonar/work/results/gipr/gipr_orthologs.csv", 0)
gipr_orthologs_fastas = mf.filter_fasta(gipr_orthologs, "/cta/users/ofkonar/work/resources/fasta/all_orthologs_pruned.fasta")
gipr_orthologs_fastas_df = mf.fasta_to_dataframe(gipr_orthologs_fastas)
gipr_consensus = mf.consensus(gipr_orthologs_fastas_df, 0.9)
gipr_consensus_posadjusted = mf.consensus_posadjusted(gipr_orthologs_fastas_df, gipr_consensus)
gipr_consensus.to_csv("/cta/users/ofkonar/work/results/gipr/gipr_consensus.csv")
gipr_consensus_posadjusted.to_csv("/cta/users/ofkonar/work/results/gipr/gipr_consensus_posadjusted.csv")

glp1r_orthologs = mf.select_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_orthologs.csv", 0)
glp1r_orthologs_fastas = mf.filter_fasta(glp1r_orthologs, "/cta/users/ofkonar/work/resources/fasta/all_orthologs_pruned.fasta")
glp1r_orthologs_fastas_df = mf.fasta_to_dataframe(glp1r_orthologs_fastas)
glp1r_consensus = mf.consensus(glp1r_orthologs_fastas_df, 0.9)
glp1r_consensus_posadjusted = mf.consensus_posadjusted(glp1r_orthologs_fastas_df, glp1r_consensus)
glp1r_consensus.to_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_consensus.csv")
glp1r_consensus_posadjusted.to_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_consensus_posadjusted.csv")

glp2r_orthologs = mf.select_csv("/cta/users/ofkonar/work/results/glp2r/glp2r_orthologs.csv", 0)
glp2r_orthologs_fastas = mf.filter_fasta(glp2r_orthologs, "/cta/users/ofkonar/work/resources/fasta/all_orthologs_pruned.fasta")
glp2r_orthologs_fastas_df = mf.fasta_to_dataframe(glp2r_orthologs_fastas)
glp2r_consensus = mf.consensus(glp2r_orthologs_fastas_df, 0.9)
glp2r_consensus_posadjusted = mf.consensus_posadjusted(glp2r_orthologs_fastas_df, glp2r_consensus)
glp2r_consensus.to_csv("/cta/users/ofkonar/work/results/glp2r/glp2r_consensus.csv")
glp2r_consensus_posadjusted.to_csv("/cta/users/ofkonar/work/results/glp2r/glp2r_consensus_posadjusted.csv")
"""
pd.set_option("display.max_rows", None, "display.max_columns", None)
"""
glr_triple_special = {26,46,48,57,93,110,112,113,138,345,385}
glr_lookout = mf.lookout("GLR", glr_triple_special)
glr_lookout.to_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_lookout.csv")

gipr_triple_special = {31,45,67,115,130,187,190,246,288,362,465}
gipr_lookout = mf.lookout("GIPR", gipr_triple_special)
gipr_lookout.to_csv("/cta/users/ofkonar/work/results/gipr/gipr_lookout.csv")

glp1r_triple_special = {31,41,52,61,70,75,92,103,105,111,123,132,137,138,140,147,156,157,160,161,171,181,195,196,202,203,204,207,208,210,212,213,216,217,224,229,230,238,242,252,258,263,276,293,294,304,305,307,311,318,327,329,338,339,342,347,362,386,389,390,397,407,413,418,430,431,432,433,434,448,452,454,456,458,459,463}
glp1r_lookout = mf.lookout("GLP1R", glp1r_triple_special)
glp1r_lookout.to_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_lookout.csv")

glp2r_triple_special = {67,109,115,117,138,192,203,215,243,252,294,299,345,360,376,381,417,423,441,454,455,526,527,534,539,545,546}
glp2r_lookout = mf.lookout("GLP2R", glp2r_triple_special)
glp2r_lookout.to_csv("/cta/users/ofkonar/work/results/glp2r/glp2r_lookout.csv")
"""

#data = pd.read_csv("/cta/users/ofkonar/work/results/all_orthologs_pruned_consensus_withhuman.csv", header = 0)
#general_consensus = data.loc[data['status'] == "C"]
#general_consensus.to_csv("/cta/users/ofkonar/work/results/general_consensus.csv")

lookfor = {143, 145, 147, 149, 151, 152, 153, 157, 163, 172, 173, 177, 179, 180, 181, 183, 184, 186, 187, 194, 231, 233, 234, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 248, 249, 251, 252, 253, 267, 270, 271, 272, 314, 315, 318, 321, 322, 325, 329, 330, 333, 343, 346, 347, 348, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 385, 391, 392, 394, 395, 396, 399, 400, 405, 406, 407, 411}
results = mf.lookout("GLR", lookfor)
print(results)
results.to_csv("/cta/users/ofkonar/work/results/wootten_significant_lookfor.csv")

