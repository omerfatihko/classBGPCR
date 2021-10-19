#read blast results from csv file
import pandas as pd

fasta_seqs = []

with open("/cta/users/ofkonar/work/fasta/gcgr_gipr_orthologs_aligned.fasta", 'r') as file:
	counter = 0
	temp = ""
	for line in file:
		if line[0] == '>':
			counter = counter + 1
			if counter == 2:
				fasta_seq = temp
				if fasta_seq != "":
					fasta_seqs.append(fasta_seq)
				temp = ""
				counter = 1				
			temp += line
		elif line[0] != ">":
			temp += line
	fasta_seq = temp
	if fasta_seq != "":
		fasta_seqs.append(fasta_seq)

human_seqs = []

for seq in fasta_seqs:
	if "HUMAN" in seq:
		human_seqs.append(seq)

list_of_fastas = []

for seq in fasta_seqs:
	split_seq = seq.split("\n")
	header = split_seq[0]
	body = "".join(map(str,split_seq[1:]))
	list_of_fastas.append({"header": header, "body": body})

list_of_humans = []

for seq in human_seqs:
	split_seq = seq.split("\n")
	header = split_seq[0]
	body = "".join(map(str,split_seq[1:]))
	list_of_humans.append({"header": header, "body": body})

doublegapindices = []

for i in range(len(list_of_humans[0]["body"])):
	if (list_of_humans[0]["body"][i] == "-") & (list_of_humans[1]["body"][i] == "-"):
		doublegapindices.append(i)

def Pruner(body, indexes):
	indexes = set(indexes)
	return "".join([char for idx, char in enumerate(body) if idx not in indexes])

for i in range(len(list_of_fastas)):
	temp = list_of_fastas[i]["body"]
	pruned_body = Pruner(temp, doublegapindices)
	list_of_fastas[i]["body"] = pruned_body


