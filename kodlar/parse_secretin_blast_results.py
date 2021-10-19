from functions import my_functions as mf
import pandas as pd

colnames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sacc", "stitle"]

#ghrhr
#read blast results
#ghrhrresults = pd.read_csv("/cta/users/ofkonar/work/results/secretin_subfamily/ghrhr/ghrhr_blast_results.csv", header = None, names = colnames)
#we get subject sequence ids (unique for every UniProt sequence, https://www.uniprot.org/help/fasta-headers)
#ghrhrids = ghrhrresults["sseqid"]
#get indices of human proteins
#indices = []
#for i, elem in enumerate(ghrhrids):
#    if 'HUMAN' in elem:
#        indices.append(i)
#target= indices[-1]
#get all protein sequences until last human sequence (last human seq will act as outgroup)
#ghrhrselected = ghrhrids[0:target+1].tolist()
#selectedghrhr = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/ghrhr/ghrhr_blast_results.csv",1)
#get fastas from unified fasta as a list
#ghrhrlist = mf.filter_fasta(ghrhrselected, "/cta/users/ofkonar/work/resources/fasta/unified_fasta.fasta")
#write the fasta list to the path
#mf.write_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/ghrhr_filtered_fastas_blast.fasta", ghrhrlist)

#pacr
#read blast results
pacrresults = pd.read_csv("/cta/users/ofkonar/work/results/secretin_subfamily/pacr/pacr_blast_results.csv", header = None, names = colnames)
#we get subject sequence ids (unique for every UniProt sequence, https://www.uniprot.org/help/fasta-headers)
pacrids = pacrresults["sseqid"]
#get indices of human proteins
pacrindices = []
for i, elem in enumerate(pacrids):
    if 'HUMAN' in elem:
        pacrindices.append(i)
pacrtarget= pacrindices[-1]
#get all protein sequences until last human sequence (last human seq will act as outgroup)
pacrselected = pacrids[0:pacrtarget+1].tolist()
#selectedpacr = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/pacr/pacr_blast_results.csv",1)
#get fastas from unified fasta as a list
pacrlist = mf.filter_fasta(pacrselected, "/cta/users/ofkonar/work/resources/fasta/unified_fasta.fasta")
#write the fasta list to the path
mf.write_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/pacr_filtered_fastas_blast.fasta", pacrlist)

#sctr
#read blast results
sctrresults = pd.read_csv("/cta/users/ofkonar/work/results/secretin_subfamily/sctr/sctr_blast_results.csv", header = None, names = colnames)
#we get subject sequence ids (unique for every UniProt sequence, https://www.uniprot.org/help/fasta-headers)
sctrids = sctrresults["sseqid"]
#get indices of human proteins
sctrindices = []
for i, elem in enumerate(sctrids):
    if 'HUMAN' in elem:
        sctrindices.append(i)
sctrtarget= sctrindices[-1]
#get all protein sequences until last human sequence (last human seq will act as outgroup)
sctrselected = sctrids[0:sctrtarget+1].tolist()
#selectedsctr = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/sctr/sctr_blast_results.csv",1)
#get fastas from unified fasta as a list
sctrlist = mf.filter_fasta(sctrselected, "/cta/users/ofkonar/work/resources/fasta/unified_fasta.fasta")
#write the fasta list to the path
mf.write_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/sctr_filtered_fastas_blast.fasta", sctrlist)

#vipr1
#read blast results
vipr1results = pd.read_csv("/cta/users/ofkonar/work/results/secretin_subfamily/vipr1/vipr1_blast_results.csv", header = None, names = colnames)
#we get subject sequence ids (unique for every UniProt sequence, https://www.uniprot.org/help/fasta-headers)
vipr1ids = vipr1results["sseqid"]
#get indices of human proteins
vipr1indices = []
for i, elem in enumerate(vipr1ids):
    if 'HUMAN' in elem:
        vipr1indices.append(i)
vipr1target= vipr1indices[-1]
#get all protein sequences until last human sequence (last human seq will act as outgroup)
vipr1selected = vipr1ids[0:vipr1target+1].tolist()
#selectedvipr1 = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/vipr1/vipr1_blast_results.csv",1)
#get fastas from unified fasta as a list
vipr1list = mf.filter_fasta(vipr1selected, "/cta/users/ofkonar/work/resources/fasta/unified_fasta.fasta")
#write the fasta list to the path
mf.write_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/vipr1_filtered_fastas_blast.fasta", vipr1list)

#vipr2
#read blast results
vipr2results = pd.read_csv("/cta/users/ofkonar/work/results/secretin_subfamily/vipr2/vipr2_blast_results.csv", header = None, names = colnames)
#we get subject sequence ids (unique for every UniProt sequence, https://www.uniprot.org/help/fasta-headers)
vipr2ids = vipr2results["sseqid"]
#get indices of human proteins
vipr2indices = []
for i, elem in enumerate(vipr2ids):
    if 'HUMAN' in elem:
        vipr2indices.append(i)
vipr2target= vipr2indices[-1]
#get all protein sequences until last human sequence (last human seq will act as outgroup)
vipr2selected = vipr2ids[0:vipr2target+1].tolist()
#selectedvipr2 = mf.select_csv("/cta/users/ofkonar/work/results/secretin_subfamily/vipr2/vipr2_blast_results.csv",1)
#get fastas from unified fasta as a list
vipr2list = mf.filter_fasta(vipr2selected, "/cta/users/ofkonar/work/resources/fasta/unified_fasta.fasta")
#write the fasta list to the path
mf.write_fasta("/cta/users/ofkonar/work/resources/fasta/secretin_subfamily/vipr2_filtered_fastas_blast.fasta", vipr2list)