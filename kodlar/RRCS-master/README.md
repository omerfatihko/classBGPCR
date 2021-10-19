# RRCS
### This code is associated with the paper from Zhou et al., "Common activation mechanism of class A GPCRs". eLife, 2019. http://dx.doi.org/10.7554/eLife.50279

RRCS (Residue-residue contact score)

RRCS can quantitatively evaluates residue contacts based on structural information. Please see workflow.jpg for the workflow of RRCS calculation.

Usage: python RRCS_calculation.py pdb_filename

Output format: Residue_1[Chain ID:Residue ID_Residue] Residue_2[Chain ID:Residue ID_Residue] RRCS

Example: for 2RH1.pdb, run "python RRCS_calculation.py 2RH1.pdb" and got output file "2RH1.pdb.cscore".

The top 10 rows of 2RH1.pdb.cscore is shown below:

A:51_ASN A:79_ASP 2.957555

A:51_ASN A:319_SER 2.449470

A:51_ASN A:80_LEU 0.746105

A:51_ASN A:323_PRO 2.886263

A:51_ASN A:76_ALA 3.645254

A:108_PHE A:170_GLN 0.400661

A:108_PHE A:166_PHE 0.976889

A:108_PHE A:112_ILE 3.025927

A:1054_THR A:1058_ILE 0.580526

A:168_PRO A:199_TYR 1.082637

