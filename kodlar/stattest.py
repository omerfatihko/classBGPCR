from functions import my_functions as mf
import pandas as pd
import time
start_time = time.time()
from os import listdir
from os.path import isfile, join
import statistics as stat
from statistics import mean
pd.set_option("display.max_rows", None) #"display.max_columns", None

csvpath = "/cta/users/ofkonar/work/results/csvs/"
#get domain info file
domaininfo = pd.read_csv("/cta/users/ofkonar/work/resources/class_B1_domains.csv", index_col = 0)
#get binary comparison files
files = [f for f in listdir(csvpath) if isfile(join(csvpath, f))]
comparisons = [f for f in files if "binary" in f]
comparisons = sorted(comparisons) #I sort the list alphabetically so it is prettier and easier to read
for file in comparisons:
	#get protein name
	prtname = file.split("_")[0]
	#get protein length from domain info
	prtlength = domaininfo.loc["C-ter", prtname.upper()]
	#keep domain lengths in one list
	domainlengths = []
	#length of N-ter (ECD)
	Nter = domaininfo.loc["N-ter", prtname.upper()]
	domainlengths.append(Nter)
	#length of 7TM
	TM7 = domaininfo.loc["TM7", prtname.upper()]
	corelength = TM7 - Nter
	domainlengths.append(corelength)
	#length of C-ter (H8 + C-ter)
	Cter = prtlength - TM7
	domainlengths.append(Cter)
	#turn the domain lengths into expected values (we will select 30 residues so normalize the length to 30)
	expected = [(x/prtlength)*30 for x in domainlengths]
	
	#chi-square
	#get the file
	comparison = pd.read_csv(csvpath + file)
	#get residues of the protein (get rid of gaps)
	withoutgaps = comparison.loc[comparison["Domain"] != "-"]
	#order the residues according to the counts (for a residue, count stands for number of binary comparisons where the residue was special(above limit frequency and different from the compared residue))
	sorteddf = withoutgaps.sort_values(by = ["count"], ascending = False)
	selected = sorteddf.head(30)
	#count the distribution of domains
	observed = []
	#get the domain column
	selecteddomains = selected["Domain"].tolist()
	sNter = 0
	s7TM = 0
	sCter = 0
	for x in selecteddomains:
		if x == "N-ter":
			sNter+=1
		elif x == "H8" or x == "C-ter":
			sCter+=1
		else:
			s7TM+=1
	observed.append(sNter)
	observed.append(s7TM)
	observed.append(sCter)
	chiresults = mf.chisquare(observed, expected)
	#print(prtname, chiresults)

	#z-score
	#keep the count of domains for samples
	countNter = []
	count7tm = []
	countCter = []
	#sample n times
	for i in range(1000):
		sample = withoutgaps.sample(n = 30)
		#count the number of domains in the list
		sampledomain = sample["Domain"].tolist()
		sampleNter = 0
		sample7TM = 0
		sampleCter = 0
		for x in sampledomain:
			if x == "N-ter":
				sampleNter+=1
			elif x == "H8" or x == "C-ter":
				sampleCter+=1
			else:
				sample7TM+=1
		countNter.append(sampleNter)
		count7tm.append(sample7TM)
		countCter.append(sampleCter)
	#calculate mean and standard deviation of domains
	meanNter = mean(countNter)
	stdNter = stat.pstdev(countNter)
	mean7tm = mean(count7tm)
	std7tm = stat.pstdev(count7tm)
	meanCter = mean(countCter)
	stdCter = stat.pstdev(countCter)
	#calculate z score for each domain type
	zsNter = (sNter - meanNter)/stdNter
	zs7tm = (s7TM - mean7tm)/std7tm
	zsCter = (sCter - meanCter)/stdCter
	#print the results
	print(prtname)
	print("Observed results: ", observed)
	print("Expected results: ", expected)
	print("sample means: ", meanNter, mean7tm, meanCter)
	print("sample st.dev.s: ", stdNter, std7tm, stdCter)
	print("Chi-square score: ", chiresults[0], ", dof: ", chiresults[1])
	print("Z-score of ECD: ", zsNter)
	print("Z-score of 7TM: ", zs7tm)
	print("Z-score of C-ter: ", zsCter)

print("My program took", time.time() - start_time, "seconds to run")