import pandas as pd
from functions import my_functions as mf
import time
start_time = time.time()

unifiedfastapath = "/cta/users/ofkonar/work/database/unified_fasta/unified_fasta.fasta"
ramporthologs = mf.select_csv("/cta/users/ofkonar/work/results/csvs/ramp3_blast_results.csv", 1)
ramporthologfastas = mf.filter_fasta(ramporthologs, unifiedfastapath)
mf.write_fasta("/cta/users/ofkonar/work/resources/ramps/ramp3_filtered_fastas_blast.fasta",ramporthologfastas)

print("My program took", time.time() - start_time, "seconds to run")