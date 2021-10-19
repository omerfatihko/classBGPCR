#read blast results from csv file
import pandas as pd

#give column names, names are arguments we gave to the blastp function
colnames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sacc", "stitle"]

#read results
results_csv = pd.read_csv("/cta/users/ofkonar/work/results/gipr/gipr_blast_results_csv.csv", header = None, names = colnames)

#we get subject sequence ids (unique for every UniProt sequence, https://www.uniprot.org/help/fasta-headers)
sseqids = results_csv["sseqid"]

#get indices of human proteins
indices = []
for i, elem in enumerate(sseqids):
    if 'HUMAN' in elem:
        indices.append(i)
#thirdind= indices[2]
fifthind= indices[4]
#print (thirdind+1) #to see how many proteins are selected
print (fifthind+1)

#get all protein sequences until 3rd human sequence (3rd human seq will act as outgroup)
#selected_seqs = sseqids[0:thirdind+1]
selected_seqs = sseqids[0:fifthind+1]

#filtering function to filter protein sequences according to results.
def FilterSeq(body, substr):
 cond = any(sub in body for sub in substr)
 if cond:
 	return body
 else:
 	return None

#open a list to append fasta sequences
fasta_seqs = []

#get fasta sequences sequence by sequence
with open('/cta/users/ofkonar/work/fasta/unified_fasta.fasta', 'r') as file:
	counter = 0
	temp = ""
	for line in file:
		if line[0] == '>':
			counter = counter + 1
			if counter == 2:
				fasta_seq = FilterSeq(temp, selected_seqs)
				if fasta_seq != None:
					fasta_seqs.append(fasta_seq)
				temp = ""
				counter = 1
			temp += line
		elif line[0] != ">":
			temp += line
	fasta_seq = FilterSeq(temp, selected_seqs)
	if fasta_seq != None:
		fasta_seqs.append(fasta_seq)

#open a file to write selected sequences
filtered_fastas=open('/cta/users/ofkonar/work/fasta/gipr_filtered_fastas_5th.fasta','w')

#write sequences to the file
for element in fasta_seqs:
     filtered_fastas.write(element)
filtered_fastas.close()