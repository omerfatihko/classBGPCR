import pandas as pd
import time
start_time = time.time()

#example script to select list of residues from a chain
"select deneme, chain A and resi 22+24+48+206"
"color white, chain A"
"hide (chain A)"
base_string = "select 11, chain R and resi "

#read the binary comparison file
comparisons = pd.read_csv("/cta/users/ofkonar/work/results/csvs/RAMPall_vs_2-3_binary_comparisons.csv")
#get residues with 14 count (different from every other GPCR)
fourteens = comparisons.loc[comparisons["count"] > 16]
#get column names
cnames = list(fourteens.columns)
#get the column name that keeps residue numbering
targetc = [x for x in cnames if "calrl_seq" in x][0]

#get the target column as a list
residues = fourteens[targetc].tolist()
print(residues)

#get the indexes of residues with + sign
indexes = ""
for ele in residues:
	onlydigits = []
	for letter in ele:
		if letter.isdigit():
			onlydigits.append(letter)
	onlydigits = "".join(onlydigits)
	onlydigits = onlydigits + "+"
	indexes = indexes + onlydigits
indexes = indexes[:-1] #last plus sign is not needed

#form the script using indexes, change the selection and chain names later, delete the plus sign at the end
script = base_string + indexes

print("My program took", time.time() - start_time, "seconds to run")

print("\n")
print(script)
print("\n")