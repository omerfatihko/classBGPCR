from functions import my_functions as mf
import pandas as pd
from os import listdir
from os.path import isfile, join
from collections import Counter


absolute_path = "/cta/users/ofkonar/work/results/"
prt_list = ["gcgr", "gipr", "glp1r", "glp2r"]

"""
for i in range(len(prt_list)):
	for j in range(i+1, len(prt_list)):
		path1 = absolute_path + prt_list[i] + "/" + prt_list[i] + "_orthologs.csv"
		path2 = absolute_path + prt_list[j] + "/" + prt_list[j] + "_orthologs.csv"
		mf.combine_csv(path1, prt_list[i], path2, prt_list[j])
"""
unified_fasta_path = "/cta/users/ofkonar/work/resources/fasta/unified_fasta.fasta"
csv_path = "/cta/users/ofkonar/work/results"
fasta_path = "/cta/users/ofkonar/work/resources/fasta/"
"""onlyfiles = [f for f in listdir(csv_path) if isfile(join(csv_path, f))]

aşağıda yaptığım şey yanlış, iterate ederken eleman mı silinir asdhfashdf

for i in onlyfiles:
	if "orthologs" not in i:
		onlyfiles.remove(i)

for i in onlyfiles:
	path = absolute_path + i
	selected = mf.select_csv(path, 0)
	fasta_list = mf.filter_fasta(selected, unified_fasta_path)
	fasta_name = i.replace(".csv", ".fasta") 
	write_path = fasta_path + fasta_name
	mf.write_fasta(write_path, fasta_list)
"""
"""
fastafiles = [f for f in listdir(fasta_path) if isfile(join(fasta_path, f))]

alignedfiles =[i for i in fastafiles if "_al" in i]

for i in alignedfiles:
	path = fasta_path + i
	truncated_i = i.replace(".fasta", "")
	fasta_file = mf.read_fasta(path)
	pruned_fasta = mf.fasta_pruner(fasta_file) 
	parts = i.split("_")
	csv_path_1 = absolute_path + parts[0] + "/" + parts[0] + "_orthologs.csv"
	csv_path_2 = absolute_path + parts[1] + "/" + parts[1] + "_orthologs.csv"
	split_pruned_dic = mf.fasta_splitter(pruned_fasta, csv_path_1, csv_path_2)
	fasta1 = split_pruned_dic["fasta_list1"]
	fasta2 = split_pruned_dic["fasta_list2"]
	write_path1 = fasta_path + parts[0] + "_" + truncated_i + "_sp.fasta"
	write_path2 = fasta_path + parts[1] + "_" + truncated_i + "_sp.fasta"
	mf.write_fasta(write_path1, fasta1)
	mf.write_fasta(write_path2, fasta2)
"""

#fastafiles = [f for f in listdir(fasta_path) if isfile(join(fasta_path, f))]
#aspfiles = [i for i in fastafiles if "_al_sp" in i]
"""
glp2r_vs_gcgr = ["/cta/users/ofkonar/work/resources/fasta/glp2r_gcgr_glp2r_orthologs_al_sp.fasta", "/cta/users/ofkonar/work/resources/fasta/gcgr_gcgr_glp2r_orthologs_al_sp.fasta"]

fasta1 = mf.read_fasta(glp2r_vs_gcgr[0])
fasta2 = mf.read_fasta(glp2r_vs_gcgr[1])

df1 = mf.fasta_to_dataframe(fasta1)
df2 = mf.fasta_to_dataframe(fasta2)

consensus1 = mf.consensus(df1, 0.9)
consensus2 = mf.consensus(df2, 0.9)

spec1 = mf.special_residues(consensus1, consensus2)
spec2 = mf.special_residues(consensus2, consensus1)

final1 = mf.finalize(df1, spec1, "GLP2R")
final2 = mf.finalize(df2, spec2, "GCGR")

final1.to_csv("/cta/users/ofkonar/work/results/glp2r/glp2r_vs_gcgr_special_residues.csv")
final2.to_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_vs_glp2r_special_residues.csv")
"""
"""
glp2r_vs_gipr = ["/cta/users/ofkonar/work/resources/fasta/glp2r_gipr_glp2r_orthologs_al_sp.fasta", "/cta/users/ofkonar/work/resources/fasta/gipr_gipr_glp2r_orthologs_al_sp.fasta"]

fasta1 = mf.read_fasta(glp2r_vs_gipr[0])
fasta2 = mf.read_fasta(glp2r_vs_gipr[1])

df1 = mf.fasta_to_dataframe(fasta1)
df2 = mf.fasta_to_dataframe(fasta2)

consensus1 = mf.consensus(df1, 0.9)
consensus2 = mf.consensus(df2, 0.9)

spec1 = mf.special_residues(consensus1, consensus2)
spec2 = mf.special_residues(consensus2, consensus1)

final1 = mf.finalize(df1, spec1, "GLP2R")
final2 = mf.finalize(df2, spec2, "GIPR")

final1.to_csv("/cta/users/ofkonar/work/results/glp2r/glp2r_vs_gipr_special_residues.csv")
final2.to_csv("/cta/users/ofkonar/work/results/gipr/gipr_vs_glp2r_special_residues.csv")
"""
"""
glp2r_vs_glp1r = ["/cta/users/ofkonar/work/resources/fasta/glp2r_glp1r_glp2r_orthologs_al_sp.fasta", "/cta/users/ofkonar/work/resources/fasta/glp1r_glp1r_glp2r_orthologs_al_sp.fasta"]

fasta1 = mf.read_fasta(glp2r_vs_glp1r[0])
fasta2 = mf.read_fasta(glp2r_vs_glp1r[1])

df1 = mf.fasta_to_dataframe(fasta1)
df2 = mf.fasta_to_dataframe(fasta2)

consensus1 = mf.consensus(df1, 0.9)
consensus2 = mf.consensus(df2, 0.9)

spec1 = mf.special_residues(consensus1, consensus2)
spec2 = mf.special_residues(consensus2, consensus1)

final1 = mf.finalize(df1, spec1, "GLP2R")
final2 = mf.finalize(df2, spec2, "GLP1R")

final1.to_csv("/cta/users/ofkonar/work/results/glp2r/glp2r_vs_glp1r_special_residues.csv")
final2.to_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_vs_glp2r_special_residues.csv")
"""
"""
glp1r_vs_gcgr = ["/cta/users/ofkonar/work/resources/fasta/glp1r_gcgr_glp1r_orthologs_al_sp.fasta", "/cta/users/ofkonar/work/resources/fasta/gcgr_gcgr_glp1r_orthologs_al_sp.fasta"]

fasta1 = mf.read_fasta(glp1r_vs_gcgr[0])
fasta2 = mf.read_fasta(glp1r_vs_gcgr[1])

df1 = mf.fasta_to_dataframe(fasta1)
df2 = mf.fasta_to_dataframe(fasta2)

consensus1 = mf.consensus(df1, 0.9)
consensus2 = mf.consensus(df2, 0.9)

spec1 = mf.special_residues(consensus1, consensus2)
spec2 = mf.special_residues(consensus2, consensus1)

final1 = mf.finalize(df1, spec1, "GLP1R")
final2 = mf.finalize(df2, spec2, "GCGR")

final1.to_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_vs_gcgr_special_residues.csv")
final2.to_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_vs_glp1r_special_residues.csv")
"""
#glp1r_vs_gipr = ["/cta/users/ofkonar/work/resources/fasta/glp1r_gipr_glp1r_orthologs_al_sp.fasta", "/cta/users/ofkonar/work/resources/fasta/gipr_gipr_glp1r_orthologs_al_sp.fasta"]

#fasta1 = mf.read_fasta(glp1r_vs_gipr[0])
#fasta2 = mf.read_fasta(glp1r_vs_gipr[1])

#df1 = mf.fasta_to_dataframe(fasta1)
#df2 = mf.fasta_to_dataframe(fasta2)

#consensus1 = mf.consensus(df1, 0.9)
#consensus1_posadjusted = mf.consensus_posadjusted(df1, consensus1)
#consensus1_posadjusted.to_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_vs_gipr_consensus.csv")
#consensus2 = mf.consensus(df2, 0.9)
#consensus2_posadjusted = mf.consensus_posadjusted(df2, consensus2)
#consensus2_posadjusted.to_csv("/cta/users/ofkonar/work/results/gipr/gipr_vs_glp1r_consensus.csv")

#spec1 = mf.special_residues(consensus1, consensus2)
#spec2 = mf.special_residues(consensus2, consensus1)

#final1 = mf.finalize(df1, spec1, "GLP1R")
#final2 = mf.finalize(df2, spec2, "GIPR")

#final1.to_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_vs_gipr_special_residues.csv")
#final2.to_csv("/cta/users/ofkonar/work/results/gipr/gipr_vs_glp1r_special_residues.csv")

for i in prt_list:
	results_path = absolute_path + i
	onlyfiles = [f for f in listdir(results_path) if isfile(join(results_path, f))]
	result_files =[f for f in onlyfiles if "_special_residues.csv" in f]
	columns =["residue"]
	for f in result_files:
		parts = f.split("_")
		columns.append(parts[2])
	columns.append("domain")
	triple_special = pd.DataFrame(columns = columns)
	df1 = pd.read_csv(results_path + "/" + result_files[0])
	df2 = pd.read_csv(results_path + "/" + result_files[1])
	df3 = pd.read_csv(results_path + "/" + result_files[2])
	res1 = df1["Name"].tolist()
	res2 = df2["Name"].tolist()
	res3 = df3["Name"].tolist()
	combined = res1 + res2 + res3
	triple = []
	for j in set(combined):
		if combined.count(j) == 3:
			triple.append(j)
	tripledf = pd.DataFrame(columns = ["res", "number"])
	for res in triple:
		tripledf = tripledf.append({"res" : res, "number" : int(''.join(x for x in res if x.isdigit()))}, ignore_index = True)
	tripledf = tripledf.sort_values('number')
	triple = tripledf["res"].tolist()
	for j in triple:
		triple_special = triple_special.append({columns[0] : j, columns[1] : df1.loc[df1["Name"] == j,"blosum80_score"].values[0], columns[2] : df2.loc[df2["Name"] == j,"blosum80_score"].values[0], columns[3] : df3.loc[df3["Name"] == j,"blosum80_score"].values[0], columns[4] : df3.loc[df3["Name"] == j,"domains"].values[0]}, ignore_index = True)
	write_path = results_path + "/" + i + "_triple_special.csv"
	triple_special.to_csv(write_path, index = False, header = True)