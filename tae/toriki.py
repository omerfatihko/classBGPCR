import time
from os import access, listdir
from os.path import isfile, join
start_time = time.time()
fastadir = "/cta/users/ofkonar/work/resources/class_B1/canonical/"
onlyfiles = [f for f in listdir(fastadir) if isfile(join(fastadir, f))]
cladefiles = [f for f in onlyfiles if "_clade.txt" in f]

with open("/cta/users/ofkonar/work/tae/raxml_parse_classB1_msa_with_outgroup.fasta") as file:
    fastalines = file.read()
    for i in cladefiles:
        prtname = i.split("_")[0].upper()
        with open(fastadir + i) as file2:
            lines = file2.readlines()
            lines = [line.rstrip() for line in lines]
            for line in lines:
                modline = prtname + "_" + line
                modline = modline.replace(" ", "_")
                modline = modline.replace(":", "")
                modline = modline.replace("(", "")
                modline = modline.replace(")", "")
                modline = modline.replace("'", "")
                fastalines = fastalines.replace(line, modline)
    fastalines = fastalines.replace(" ", "_")
    fastalines = fastalines.replace(":", "")
    fastalines = fastalines.replace("(", "")
    fastalines = fastalines.replace(")", "")
    fastalines = fastalines.replace("'", "")
with open("/cta/users/ofkonar/work/tae/raxml_parsable.fasta", "w") as file3:
    file3.write(fastalines)
    

print("My program took", time.time() - start_time, "seconds to run")