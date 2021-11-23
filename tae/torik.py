import time
from os import access, listdir
from os.path import isfile, join
start_time = time.time()
fastadir = "/cta/users/ofkonar/work/resources/class_B1/canonical/"
onlyfiles = [f for f in listdir(fastadir) if isfile(join(fastadir, f))]
cladefiles = [f for f in onlyfiles if "_clade.txt" in f]

with open("/cta/users/ofkonar/work/tae/fasttreestr.txt") as file:
    newickline = file.read()
    for i in cladefiles:
        prtname = i.split("_")[0].upper()
        with open(fastadir + i) as file2:
            lines = file2.readlines()
            lines = [line.rstrip() for line in lines]
            for line in lines:
                modline = prtname + "_" + line
                newickline = newickline.replace(line, modline)
with open("/cta/users/ofkonar/work/tae/fasttreestrnew.txt", "w") as file3:
    file3.write(newickline)
        
print("My program took", time.time() - start_time, "seconds to run")