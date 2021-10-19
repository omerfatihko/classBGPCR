import pandas as pd

pd.set_option("display.max_rows", None, "display.max_columns", None)

a = pd.read_csv("/cta/users/ofkonar/work/kodlar/RRCS-master/structures/5yqz_glucagon_active.pdb.cscore", delim_whitespace=True, header=None, engine="python")

dims = a.shape
gcgr_inactive = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"])

for i in range(dims[0]):
	row = a.iloc[i]
	chain_1 = row[0].split(":")[0]
	chain_2 = row[1].split(":")[0]
	if chain_1 == "R" and chain_2 == "R":
		gcgr_inactive = gcgr_inactive.append({"chain_1" : chain_1, "residue_1" : row[0].split(":")[1], "chain_2" : chain_2, "residue_2" : row[1].split(":")[1], "RRCS" : row[2]}, ignore_index = True)

b = pd.read_csv("/cta/users/ofkonar/work/kodlar/RRCS-master/structures/6lmk.pdb.cscore", delim_whitespace = True, engine = "python" , header=None)

dims2 = b.shape
gcgr_active = pd.DataFrame(columns = ["chain_1", "residue_1", "chain_2", "residue_2", "RRCS"])

for i in range(dims2[0]):
	row = b.iloc[i]
	chain_1 = row[0].split(":")[0]
	chain_2 = row[1].split(":")[0]
	if chain_1 == "R" and chain_2 == "R":
		gcgr_active = gcgr_active.append({"chain_1" : chain_1, "residue_1" : row[0].split(":")[1], "chain_2" : chain_2, "residue_2" : row[1].split(":")[1], "RRCS" : row[2]}, ignore_index = True)

pinactive = gcgr_inactive[["residue_1", "residue_2"]].to_records(index = False).tolist()
inactive = []
for ele in pinactive:
	a = int(''.join(x for x in ele[0] if x.isdigit()))
	b = int(''.join(x for x in ele[1] if x.isdigit()))
	if a <478 and b<478:
		if a < b:
			inactive.append(ele)
		else:
			inactive.append(reversed(ele))

pactive = gcgr_active[["residue_1", "residue_2"]].to_records(index = False).tolist()
active = []
for ele in pactive:
	a = int(''.join(x for x in ele[0] if x.isdigit()))
	b = int(''.join(x for x in ele[1] if x.isdigit()))
	if a <478 and b<478:
		if a < b:
			active.append(ele)
		else:
			active.append(reversed(ele))

switch_off = []
switch_on = []
repack = []
for inactive_ele in inactive:
 if inactive_ele in active:
 	repack.append(inactive_ele)
 else:
 	switch_off.append(inactive_ele)

for active_ele in active:
	if active_ele not in inactive:
		switch_on.append(active_ele)

switch_off_df = pd.DataFrame(columns = ["residue_1", "residue_2", "dRRCS"])
for ele in switch_off:
	residue_1 = ele[0]
	residue_2 = ele[1]
	dRRCS = gcgr_inactive["RRCS"].loc[(gcgr_inactive["residue_1"] == residue_1) & (gcgr_inactive["residue_2"] == residue_2)].values[0]
	switch_off_df = switch_off_df.append({"residue_1" : residue_1, "residue_2" : residue_2, "dRRCS" : dRRCS}, ignore_index = True)

switch_on_df = pd.DataFrame(columns = ["residue_1", "residue_2", "dRRCS"])
for ele in switch_on:
	residue_1 = ele[0]
	residue_2 = ele[1]
	dRRCS = gcgr_active["RRCS"].loc[(gcgr_active["residue_1"] == residue_1) & (gcgr_active["residue_2"] == residue_2)].values[0]
	switch_on_df = switch_on_df.append({"residue_1" : residue_1, "residue_2" : residue_2, "dRRCS" : dRRCS}, ignore_index = True)

repack_df = pd.DataFrame(columns = ["residue_1", "residue_2", "dRRCS"])
for ele in repack:
	residue_1 = ele[0]
	residue_2 = ele[1]
	RRCS1 = gcgr_inactive["RRCS"].loc[(gcgr_inactive["residue_1"] == residue_1) & (gcgr_inactive["residue_2"] == residue_2)].values[0]
	RRCS2 = gcgr_active["RRCS"].loc[(gcgr_active["residue_1"] == residue_1) & (gcgr_active["residue_2"] == residue_2)].values[0]
	repack_df = repack_df.append({"residue_1" : residue_1, "residue_2" : residue_2, "dRRCS" : RRCS2 - RRCS1}, ignore_index = True)


