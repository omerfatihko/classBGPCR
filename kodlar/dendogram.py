from functions import my_functions as mf
import pandas as pd

#pd.set_option("display.max_rows", None, "display.max_columns", None)

clade = []
with open("/cta/users/ofkonar/work/results/glp2r/glp2r_clade.txt", 'r') as file:
	for line in file:
		clade.append(line)

a = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/glp2r_glp1r_glp2r_orthologs_al_sp.fasta")
adf = mf.fasta_to_dataframe(a)

special = pd.read_csv("/cta/users/ofkonar/work/results/glp2r/glp2r_triple_special.csv", header = 0)

new_newick = mf.dendogram_seq(adf, special, clade[2])
clade[2] = new_newick

new_txt = open('/cta/users/ofkonar/work/tae/glp2r_newick_with_special.txt','w')

for element in clade:
     new_txt.write(element)
new_txt.close()