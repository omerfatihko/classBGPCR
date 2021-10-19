#OBSOLETE CODES

def double_consensus(lst1, lst2):
	"""This function takes 2 equal length consensus lists and make a bigger consensus list.
	For double consensus list two consensus lists should agree at every position otherwise
	that position is taken as gap"""

	double_consensus = []
	for i in range(len(lst1)):
		if lst1[i] == lst2[i]:
			double_consensus.append(lst1[i])
		else:
			double_consensus.append(["-"])

	return double_consensus

def pos_adjust(df, lst):
	"""This function takes a list of residues and adjust their position to the original
	human sequence"""

	cnames = list(df.columns)
	cname = ""
	for element in cnames:
		if element.find("HUMAN") != -1:
			cname = element	
	human_seq = df[cname].tolist()
	indexes = []
	for i in range(len(human_seq)):
		if human_seq[i] == "-":
			indexes.append(i)

	for ele in sorted(indexes, reverse = True):  
		del human_seq[ele]
	for ele in sorted(indexes, reverse = True):  
		del lst[ele]

	residue_df = pd.DataFrame(list(zip(human_seq, lst)), columns =['human', 'Spec_res'])

	return residue_df

def special_res(df):
	"""This function takes the resulting dataframeof pos_adjuct() function and returns list of 
	special residues"""

	special_res = []
	dims = df.shape[0]
	for i in range(dims):
		if df.iloc[i,1][0] != "-":
			special_res.append(df.iloc[i,1][0]+str(i+1))

	return special_res

def res_to_Wootten(special_res, prt_name):
	"""This function translates residue numbers to Wootten numbering scheme
	Column names : 'CALCR', 'CALCRL', 'CRHR1', 'CRHR2', 'GHRHR', 'GIPR',
	'GLP1R', 'GLP2R', 'GCGR', 'SCTR', 'PTH1R', 'PTH2R', 'PAC1R','VPAC1R', 'VPAC2R'"""

	Wootten = []
	df = pd.read_csv("/cta/users/ofkonar/work/resources/residue_table.csv")
	for res in special_res:
		row = df.loc[df[prt_name] == res]
		if row.empty == False:
			Wootten.append(row["Wootten"].values[0])
		else:
			Wootten.append("-")

	res_to_Wootten = pd.DataFrame(list(zip(special_res, Wootten)), columns =['Residues', 'Wootten'])

	return res_to_Wootten 

####BUNA DEVAM ETMEYE ÇALIŞALIM
string = "((tr|A0A1S3NQC3|A0A1S3NQC3_SALSA:0.0208508846,tr|A0A060WGN5|A0A060WGN5_ONCMY:0.0642213927)79.7/100:0.0119732114,(tr|A0A1S3SZ13|A0A1S3SZ13_SALSA:0.0216563050,tr|A0A060W6Z9|A0A060W6Z9_ONCMY:0.0290970039)74.2/100:0.0131374042)99.9/100:0.1401977352"

file_to_read = open("/cta/users/ofkonar/work/results/gipr/gipr_iqtree_2/gipr_iqtree_boostrap_modelfinder_2.treefile", "r")
contents = file_to_read.read()
file_to_read.close

istart = []  # stack of indices of opening parentheses
d = {}

for i, c in enumerate(contents):
    if c == '(':
         istart.append(i)
    if c == ')':
        try:
            d[istart.pop()] = i
        except IndexError:
            print('Too many closing parentheses')
if istart:  # check if stack is empty afterwards
    print('Too many opening parentheses')
print(d)