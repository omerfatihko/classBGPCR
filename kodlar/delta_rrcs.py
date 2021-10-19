import pandas as pd

gcgr_active = pd.read_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_RRCS_active.csv", header = 0, index_col = 0)
gcgr_inactive = pd.read_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_RRCS_inactive.csv", header = 0, index_col = 0)
glp1r_active = pd.read_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_RRCS_active.csv", header = 0, index_col = 0)
glp1r_inactive = pd.read_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_RRCS_inactive.csv", header = 0, index_col = 0)

gcgr_active_columns = list(gcgr_active.columns)
gcgr_inactive_columns = list(gcgr_inactive.columns)
glp1r_active_columns = list(glp1r_active.columns)
glp1r_inactive_columns = list(glp1r_inactive.columns)

gcgr_delta_columns = []
for col1 in gcgr_active_columns:
	for col2 in gcgr_inactive_columns:
		newcol = col1 + "vs" + col2
		gcgr_delta_columns.append(newcol)

glp1r_delta_columns = []
for col1 in glp1r_active_columns:
	for col2 in glp1r_inactive_columns:
		newcol = col1 + "vs" + col2
		glp1r_delta_columns.append(newcol)

gcgr_delta_indexes = list(gcgr_active.index)
glp1r_delta_indexes = list(glp1r_active.index)

gcgr_delta = pd.DataFrame(columns = gcgr_delta_columns, index = gcgr_delta_indexes)
glp1r_delta = pd.DataFrame(columns = glp1r_delta_columns, index = glp1r_delta_indexes)

for cols in gcgr_delta.columns:
	part1 = cols.split("vs")[0]
	part2 = cols.split("vs")[1]
	gcgr_delta[cols] = gcgr_active[part1] - gcgr_inactive[part2]

for cols in glp1r_delta.columns:
	part1 = cols.split("vs")[0]
	part2 = cols.split("vs")[1]
	glp1r_delta[cols] = glp1r_active[part1] - glp1r_inactive[part2]

gcgr_switch_on = []
gcgr_switch_off = []
gcgr_repack = []
gcgr_dud = []
for i in range(gcgr_delta.shape[0]): 
	active_list = gcgr_active.iloc[i,:].values.tolist()
	inactive_list = gcgr_inactive.iloc[i,:].values.tolist()
	switchon_count = 0
	switchoff_count = 0
	repack_count = 0
	dud_count = 0
	for j in active_list:
		for k in inactive_list:
			if j == 0 and k == 0:
				dud_count += 1
			elif j == 0 and k != 0:
				switchoff_count += 1
			elif j != 0 and k == 0:
				switchon_count += 1
			else:
				repack_count +=1
	gcgr_switch_on.append(switchon_count)
	gcgr_switch_off.append(switchoff_count)
	gcgr_repack.append(repack_count)
	gcgr_dud.append(dud_count)

gcgr_delta["switch_off"] = gcgr_switch_off
gcgr_delta["switch_on"] = gcgr_switch_on
gcgr_delta["repack"] = gcgr_repack
gcgr_delta["dud"] = gcgr_dud

glp1r_switch_on = []
glp1r_switch_off = []
glp1r_repack = []
glp1r_dud = []
for i in range(glp1r_delta.shape[0]): 
	active_list = glp1r_active.iloc[i,:].values.tolist()
	inactive_list = glp1r_inactive.iloc[i,:].values.tolist()
	switchon_count = 0
	switchoff_count = 0
	repack_count = 0
	dud_count = 0
	for j in active_list:
		for k in inactive_list:
			if j == 0 and k == 0:
				dud_count += 1
			elif j == 0 and k != 0:
				switchoff_count += 1
			elif j != 0 and k == 0:
				switchon_count += 1
			else:
				repack_count +=1
	glp1r_switch_on.append(switchon_count)
	glp1r_switch_off.append(switchoff_count)
	glp1r_repack.append(repack_count)
	glp1r_dud.append(dud_count)

glp1r_delta["switch_off"] = glp1r_switch_off
glp1r_delta["switch_on"] = glp1r_switch_on
glp1r_delta["repack"] = glp1r_repack
glp1r_delta["dud"] = glp1r_dud

#gcgr_delta.to_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_RRCS_delta.csv")
#glp1r_delta.to_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_RRCS_delta.csv")

gcgr_delta_critic = gcgr_delta
glp1r_delta_critic = glp1r_delta

indexes_to_drop = list(gcgr_delta_critic.index)
idx = []
for i in range(gcgr_delta_critic.shape[0]):
	row = gcgr_delta_critic.iloc[i,:].values.tolist()
	a = row[:len(row)-4]
	if all(val > 0.2 for val in a):
		idx.append(i)
	elif all(val < -0.2 for val in a):
		idx.append(i)

for ele in sorted(idx, reverse = True):
	del indexes_to_drop[ele]

gcgr_delta_critic = gcgr_delta_critic.drop(indexes_to_drop)
#gcgr_delta_critic.to_csv("/cta/users/ofkonar/work/results/gcgr/gcgr_RRCS_delta_critic.csv")

indexes_to_drop = list(glp1r_delta_critic.index)
idx = []
for i in range(glp1r_delta_critic.shape[0]):
	row = glp1r_delta_critic.iloc[i,:].values.tolist()
	a = row[:len(row)-4]
	if all(val > 0.2 for val in a):
		idx.append(i)
	elif all(val < -0.2 for val in a):
		idx.append(i)

for ele in sorted(idx, reverse = True):
	del indexes_to_drop[ele]

glp1r_delta_critic = glp1r_delta_critic.drop(indexes_to_drop)
#glp1r_delta_critic.to_csv("/cta/users/ofkonar/work/results/glp1r/glp1r_RRCS_delta_critic.csv")