from functions import my_functions as mf
import pandas as pd

pd.set_option("display.max_rows", None, "display.max_columns", None)

a = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/older_fastas/gcgr_split_pruned.fasta")
b = mf.read_fasta("/cta/users/ofkonar/work/resources/fasta/older_fastas/gipr_split_pruned.fasta")

adf = mf.fasta_to_dataframe(a)
bdf = mf.fasta_to_dataframe(b)

acons = mf.consensus(adf, 0.9)
#acons_posadjusted = mf.consensus_posadjusted(adf, acons)
#acons_posadjusted.to_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_consensus.csv")

bcons = mf.consensus(bdf, 0.9)
#bcons_posadjusted = mf.consensus_posadjusted(bdf, bcons)
#bcons_posadjusted.to_csv("/cta/users/ofkonar/work/results/gipr/gipr_consensus.csv")

aspec = mf.special_residues(acons, bcons)
#afinal = mf.finalize(adf, aspec, "GCGR")
#afinal_ordered = afinal.sort_values(by = ["blosum80_score", "frequency"], ascending = [True, False])
#afinal_ordered.to_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_special_residues_ordered.csv")

#for GIPR
#bspec = mf.special_residues(bcons, acons)
#bfinal = mf.finalize(bdf, bspec, "GIPR")
#bfinal_ordered = bfinal.sort_values(by = ["blosum80_score", "frequency"], ascending = [True, False])
#bfinal_ordered.to_csv("/cta/users/ofkonar/work/results/gipr/gipr_special_residues_ordered.csv")

gapless = mf.dendogram_seq(adf)
print(gapless)