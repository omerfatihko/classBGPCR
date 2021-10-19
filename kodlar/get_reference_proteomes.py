import urllib.request
import csv
import time
start_time = time.time()
#Proteome ID file accesssion date and url
"14.07.2021 https://www.uniprot.org/proteomes/?query=*&fil=taxonomy%3A%22Eukaryota+%5B2759%5D%22+AND+reference%3Ayes"

#provide raw URL, XXX will be replaced with proteome ID
raw_url = 'https://www.uniprot.org/uniprot/?query=proteome:XXX&format=fasta'

#read the tab delimited file containing proteome IDs and organism IDs
with open("/cta/users/ofkonar/work/database/14._07_2021_proteomes-filtered-taxonomy__Eukaryota+[2759]_+AND+reference_yes.tab", newline="") as file:
    file_reader = csv.DictReader(file, delimiter = "\t")
    for line in file_reader:
    	#get required identifiers
        proteome_id = line["Proteome ID"]
        organism = line["Organism"].replace(" ", "_") #remove gaps
        organism = organism.replace("/", "-") #/ sign is considered as directory
        organism = organism[:24] #too long names cause errors
        organism_id = line["Organism ID"]
        #provide url to access reference proteome
        currenturl = raw_url.replace("XXX", proteome_id)
        #provide fasta file name (includes important identifiers)
        fastaname = proteome_id + "_" + organism_id + "_" + organism + ".fasta"
        #provide path to write down fastas
        writepath = "/cta/users/ofkonar/work/database/reference_proteomes/"
        #request the proteome, write it down to the provided path
        urllib.request.urlretrieve(currenturl,writepath + fastaname)

print("My program took", time.time() - start_time, "seconds to run")