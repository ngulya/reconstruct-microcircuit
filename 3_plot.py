import matplotlib as mpl
import matplotlib.pyplot as plt
import json
import os
import numpy as np
import time
import shutil
import random
import seaborn as sns
#L5->L6->L4->L23->L6
def bursting_per_time(layer_nums, bursting, folder=''):
	# layer_nums = {'L1': 46, 'L23': 612, 'L4': 704, 'L5': 908, 'L6': 1891}
	##########bad syn
	# bursting = {160: {'L6': 29, 'L5': 15, 'L4': 8, 'L23': 5, 'L1': 2}, 170: {'L5': 10, 'L6': 7, 'L4': 5, 'L23': 3, 'L1': 2}, 180: {'L6': 27, 'L5': 15, 'L1': 5, 'L4': 4, 'L23': 2}, 190: {'L6': 20, 'L5': 8, 'L1': 5, 'L23': 3, 'L4': 2}, 200: {'L5': 11, 'L6': 10, 'L4': 7, 'L1': 3, 'L23': 2}, 210: {'L5': 7, 'L6': 6, 'L23': 3, 'L4': 2, 'L1': 0}, 220: {'L5': 17, 'L6': 11, 'L4': 5, 'L1': 4, 'L23': 3}, 230: {'L23': 82, 'L6': 10, 'L5': 7, 'L4': 3, 'L1': 1}, 240: {'L4': 636, 'L23': 250, 'L5': 57, 'L6': 26, 'L1': 11}, 250: {'L6': 1893, 'L5': 792, 'L4': 96, 'L23': 16, 'L1': 3}, 260: {'L5': 152, 'L4': 31, 'L6': 30, 'L1': 3, 'L23': 2}, 270: {'L5': 201, 'L6': 90, 'L4': 35, 'L23': 7, 'L1': 0}, 280: {'L6': 585, 'L5': 248, 'L4': 13, 'L23': 7, 'L1': 1}, 290: {'L6': 1251, 'L4': 155, 'L5': 105, 'L23': 7, 'L1': 2}, 300: {'L6': 174, 'L5': 136, 'L4': 124, 'L23': 5, 'L1': 0}, 310: {'L4': 69, 'L6': 24, 'L5': 12, 'L23': 3, 'L1': 1}, 320: {'L6': 33, 'L5': 23, 'L4': 22, 'L23': 4, 'L1': 1}, 330: {'L6': 32, 'L5': 14, 'L4': 11, 'L23': 5, 'L1': 0}, 340: {'L5': 22, 'L6': 12, 'L4': 6, 'L1': 4, 'L23': 4}, 350: {'L5': 12, 'L6': 7, 'L4': 6, 'L23': 5, 'L1': 3}}
	# bursting = {1320: {'L5': 7, 'L4': 5, 'L6': 2, 'L1': 1, 'L23': 1}, 1330: {'L6': 11, 'L5': 10, 'L23': 4, 'L4': 4, 'L1': 1}, 1340: {'L5': 13, 'L6': 10, 'L23': 3, 'L4': 2, 'L1': 0}, 1350: {'L5': 9, 'L4': 8, 'L6': 6, 'L23': 5, 'L1': 1}, 1360: {'L6': 15, 'L5': 12, 'L4': 6, 'L23': 4, 'L1': 1}, 1370: {'L5': 13, 'L23': 3, 'L4': 3, 'L6': 3, 'L1': 1}, 1380: {'L5': 11, 'L6': 9, 'L4': 4, 'L23': 2, 'L1': 1}, 1390: {'L5': 14, 'L6': 13, 'L4': 4, 'L23': 3, 'L1': 0}, 1400: {'L5': 30, 'L6': 19, 'L23': 6, 'L4': 3, 'L1': 1}, 1410: {'L6': 113, 'L5': 34, 'L4': 8, 'L23': 3, 'L1': 1}, 1420: {'L6': 1691, 'L5': 337, 'L4': 9, 'L23': 4, 'L1': 0}, 1430: {'L6': 1214, 'L4': 662, 'L5': 532, 'L23': 248, 'L1': 1}, 1440: {'L5': 653, 'L6': 320, 'L23': 211, 'L4': 106, 'L1': 4}, 1450: {'L6': 82, 'L5': 50, 'L4': 29, 'L23': 3, 'L1': 1}, 1460: {'L6': 27, 'L5': 17, 'L4': 10, 'L23': 4, 'L1': 1}, 1470: {'L5': 36, 'L6': 20, 'L23': 10, 'L4': 8, 'L1': 1}, 1480: {'L5': 18, 'L6': 14, 'L4': 5, 'L23': 3, 'L1': 2}, 1490: {'L6': 16, 'L5': 15, 'L4': 10, 'L23': 5, 'L1': 1}, 1500: {'L5': 12, 'L6': 12, 'L4': 6, 'L23': 4, 'L1': 0}, 1510: {'L5': 16, 'L6': 11, 'L23': 6, 'L4': 6, 'L1': 0}}
	# bursting = {2450: {'L5': 12, 'L6': 8, 'L4': 6, 'L23': 4, 'L1': 0}, 2460: {'L5': 18, 'L6': 7, 'L23': 5, 'L4': 3, 'L1': 0}, 2470: {'L5': 12, 'L6': 8, 'L23': 2, 'L1': 1, 'L4': 1}, 2480: {'L5': 19, 'L4': 6, 'L6': 4, 'L23': 3, 'L1': 1}, 2490: {'L5': 23, 'L6': 12, 'L4': 9, 'L23': 7, 'L1': 0}, 2500: {'L5': 22, 'L6': 19, 'L23': 5, 'L4': 4, 'L1': 0}, 2510: {'L6': 86, 'L5': 30, 'L23': 2, 'L4': 1, 'L1': 0}, 2520: {'L6': 1661, 'L5': 317, 'L4': 13, 'L23': 9, 'L1': 0}, 2530: {'L6': 1012, 'L4': 661, 'L5': 604, 'L23': 204, 'L1': 2}, 2540: {'L5': 569, 'L6': 387, 'L23': 193, 'L4': 100, 'L1': 2}, 2550: {'L6': 212, 'L5': 56, 'L4': 27, 'L23': 6, 'L1': 3}, 2560: {'L6': 33, 'L5': 29, 'L4': 16, 'L23': 4, 'L1': 1}, 2570: {'L5': 28, 'L6': 20, 'L4': 7, 'L23': 5, 'L1': 0}, 2580: {'L5': 18, 'L6': 12, 'L23': 7, 'L4': 6, 'L1': 1}, 2590: {'L5': 11, 'L6': 10, 'L4': 8, 'L23': 5, 'L1': 0}, 2600: {'L5': 17, 'L6': 15, 'L23': 5, 'L4': 5, 'L1': 1}, 2610: {'L5': 9, 'L6': 8, 'L4': 6, 'L23': 3, 'L1': 0}, 2620: {'L6': 8, 'L4': 7, 'L5': 7, 'L23': 6, 'L1': 1}, 2630: {'L5': 10, 'L4': 8, 'L6': 6, 'L23': 5, 'L1': 0}, 2640: {'L5': 8, 'L6': 4, 'L23': 3, 'L1': 0, 'L4': 0}}
	
	##########norm syn
	# bursting = {210: {'L5': 8, 'L6': 8, 'L1': 2, 'L4': 2, 'L23': 1}, 220: {'L5': 16, 'L23': 11, 'L6': 7, 'L4': 4, 'L1': 0}, 230: {'L23': 199, 'L4': 18, 'L6': 13, 'L5': 11, 'L1': 6}, 240: {'L6': 1217, 'L4': 649, 'L5': 636, 'L23': 142, 'L1': 9}, 250: {'L6': 773, 'L5': 227, 'L4': 83, 'L23': 7, 'L1': 2}, 260: {'L5': 161, 'L6': 47, 'L4': 24, 'L23': 5, 'L1': 0}, 270: {'L5': 256, 'L6': 174, 'L4': 35, 'L23': 8, 'L1': 5}, 280: {'L6': 860, 'L5': 216, 'L4': 42, 'L23': 8, 'L1': 1}, 290: {'L6': 960, 'L4': 155, 'L5': 142, 'L23': 7, 'L1': 0}, 300: {'L4': 116, 'L5': 48, 'L6': 47, 'L1': 2, 'L23': 0}}
	# bursting = {1240: {'L6': 25, 'L5': 18, 'L4': 7, 'L23': 3, 'L1': 0}, 1250: {'L6': 48, 'L5': 19, 'L4': 6, 'L1': 0, 'L23': 0}, 1260: {'L6': 155, 'L5': 37, 'L23': 5, 'L4': 4, 'L1': 0}, 1270: {'L6': 1648, 'L5': 362, 'L4': 21, 'L23': 6, 'L1': 2}, 1280: {'L6': 1005, 'L4': 644, 'L5': 542, 'L23': 26, 'L1': 4}, 1290: {'L5': 422, 'L23': 415, 'L6': 193, 'L4': 72, 'L1': 2}, 1300: {'L5': 206, 'L6': 180, 'L4': 73, 'L23': 13, 'L1': 5}, 1310: {'L5': 73, 'L6': 48, 'L4': 7, 'L23': 3, 'L1': 1}, 1320: {'L6': 21, 'L5': 19, 'L4': 9, 'L23': 7, 'L1': 0}, 1330: {'L6': 23, 'L5': 16, 'L4': 10, 'L23': 7, 'L1': 0}}
	# bursting = {2220: {'L5': 17, 'L6': 13, 'L23': 4, 'L4': 3, 'L1': 1}, 2230: {'L5': 24, 'L6': 12, 'L4': 5, 'L23': 3, 'L1': 0}, 2240: {'L5': 36, 'L6': 32, 'L23': 6, 'L4': 2, 'L1': 1}, 2250: {'L6': 737, 'L5': 140, 'L4': 10, 'L23': 4, 'L1': 0}, 2260: {'L6': 1497, 'L5': 730, 'L4': 557, 'L23': 12, 'L1': 3}, 2270: {'L23': 429, 'L6': 421, 'L5': 338, 'L4': 123, 'L1': 0}, 2280: {'L5': 306, 'L6': 306, 'L4': 109, 'L23': 67, 'L1': 4}, 2290: {'L6': 115, 'L5': 49, 'L4': 20, 'L23': 2, 'L1': 0}, 2300: {'L6': 34, 'L5': 33, 'L4': 12, 'L23': 8, 'L1': 2}, 2310: {'L5': 23, 'L6': 22, 'L4': 20, 'L23': 6, 'L1': 1}}
	l = len(bursting)
	plt.figure(figsize=(7,14))
	# plt.figure(figsize=(5,8))
	mpl.rcParams['font.size'] = 13.0
	for i, (time, dict_num_spikes) in enumerate(bursting.items()):
		ax = plt.subplot(l, 1, i+1, frameon=False)
		relative_nums = {}
		for layer in layer_nums:
			relative_nums[layer] = dict_num_spikes[layer]/layer_nums[layer]	

		sum_alls = sum(relative_nums.values())
		for layer in relative_nums:
			relative_nums[layer] = relative_nums[layer]/sum_alls

		barlist = ax.bar(list(relative_nums.keys()), list(relative_nums.values()), color='blue')
		
		lnums = list(relative_nums.values())
		barlist[lnums.index(max(lnums))].set_color('red')
		# for i, (k, v) in enumerate(relative_nums.items()):
		# 	ax.text(k, 0, str(np.round(v, 3)), color='green')
		# for idx,rect in enumerate(barlist):
		# 	height = rect.get_height()
		# 	height = 0.05*height
		# 	ax.text(rect.get_x() + rect.get_width()/2., height, np.round(list(relative_nums.values())[idx], 3), ha='center', va='bottom', color = 'blue', rotation=0)

		ax.set_yticklabels([time])
		ax.set_ylim([0, 1])
		if i < l - 1:
			ax.set_xticklabels([''])
	if folder == '':
		plt.show()
	else:
		plt.savefig(opj(folder, 'bursting_per_time.png'))
	# exit('-')
# bursting_per_time()

def get_neuron_number(folder='template'):
	num_neurons = 0
	for _template_paths in os.listdir(folder):
		if _template_paths.endswith('.hoc'):
			num_neurons += 1
	return num_neurons

def plot_all_cells_parallel(folder='template'):
	num_cells = get_neuron_number(folder)
	t = []
	l_soma_v = []
	l_name = []


	for file in os.listdir('result'):
		if file.startswith('host_') and file.startswith('host_') and file.endswith('.json'):
			with open(f'result/{file}', 'r') as outfile:
				result = json.load(outfile)
			if not len(t):
				t = result['t']
			for cell_id, soma_v in result['record_soma_v'].items():
				l_name.append(cell_id)
				l_soma_v.append(soma_v)

	for i in range(num_cells):
		label = l_name[i]
		if i == 3:
			label = 'Main: ' + label
		plt.plot(t, l_soma_v[i], label = label)

	plt.legend()
	plt.show()

def get_inh_exi():
	global should_mtypes
	inhib = {}
	excit = {}
	for t in should_mtypes:
		inhib[t] = 0
		excit[t] = 0
	
	with open('BBPjson/pathways_physiology_factsheets_simplified.json', 'r') as fd:
		pathways = json.load(fd)

	for k, v in pathways.items():
		cell = k.split(':')[0].replace('-', '_')
		if v['synapse_type'].find('Inhibitory') >= 0:
			inhib[cell] += 1
		else:
			excit[cell] += 1

	for k in should_mtypes:
		if inhib[k] > excit[k]:
			inhib[k] = 1
			excit[k] = 0
		else:
			inhib[k] = 0
			excit[k] = 1

	for k in should_mtypes:
		if inhib[k] == 1:
			print(k, 'inh')
		else:
			print(k, 'exc')
	return inhib, excit


def get_sorted_cellids(cell_ids, should_mtypes):
	wascellid = len(cell_ids)
	dict_sorted_cell = {}
	for i in should_mtypes:
		dict_sorted_cell[i] = []
		for cell in cell_ids:
			if cell.startswith(i):
				dict_sorted_cell[i].append(cell)
		for cell in dict_sorted_cell[i]:
			cell_ids.remove(cell)
	sorted_cell = []
	for k, v in dict_sorted_cell.items():
		sorted_cell.extend(v)
	if len(cell_ids) != 0:
		exit('-155')
	return sorted_cell, wascellid

def create_fig30(inhib, excit, res_folder):
	global should_mtypes, should_etypes

	if not os.path.exists(opj(res_folder, 'cell.json')):
		return False

	print(opj(res_folder, 'cell.json'))
	with open(opj(res_folder, 'cell.json'), 'r') as fd:
		cell_ids = json.load(fd)['cell_ids']

	sq_ex = [[], ['PC'], ['PC', 'SS', 'SP'], ['STPC', 'UTPC', 'TTPC2', 'TTPC1'], ['TPC_L4', 'BPC', 'IPC', 'TPC_L1', 'UTPC', 'TPC_L4']]
	sq_in = [['NGC_SA', 'NGC_DA', 'SLAC', 'HAC', 'DLAC', 'DAC'], ['BP', 'SBC', 'NBC', 'LBC', 'NGC', 'MC', 'DBC', 'ChC', 'BTC'], ['BP', 'SBC', 'NBC', 'LBC', 'NGC', 'MC', 'DBC', 'ChC', 'BTC'], ['BP', 'SBC', 'NBC', 'LBC', 'NGC', 'MC', 'DBC', 'ChC', 'BTC'], ['BP', 'SBC', 'NBC', 'LBC', 'DBC', 'NGC', 'MC', 'ChC', 'BTC']]
	# plt.figure(figsize=(16,9))
	angelex = [0, -37, -37, 110, 130]
	angelin = [100, 100, 100, 80, 80]
	
	mpl.rcParams['font.size'] = 13.0
	# plt.figure(figsize=(6, 8))
	plt.figure(figsize=(18, 16))
	all_heat = []
	for _i, layer in enumerate(['L1', 'L23', 'L4', 'L5', 'L6']):
		i, e = 0, 0
		print('\n', layer)
		mtype_inhib = {}
		mtype_excit = {}
		numl = 0
		heatlist = []
		this_layer_mtypes = []

		for cell in cell_ids:
			cl = cell.split('_')[0]
			if cl != layer:
				continue
			for mtype in should_mtypes:
				if cell.startswith(mtype):
					break
			for etype in should_etypes:
				if etype in cell:
					break
			# print(cell, mtype, etype)
			i += inhib[mtype]
			e += excit[mtype]
			withoutlay = mtype.replace(layer+'_', '')
			if inhib[mtype] == 1:
				mtype_inhib.setdefault(withoutlay, 0)
				mtype_inhib[withoutlay] += 1
			else:
				mtype_excit.setdefault(withoutlay, 0)
				mtype_excit[withoutlay] += 1

		print(layer, f' |  i:{np.round(i/(i+e), 3)} e:{np.round(e/(i+e), 3)}	| i:{i}  e:{e}')
		status, newmtype_inhib = dict_sort_by_seq(mtype_inhib, sq_in[_i])
		if not status: exit('er1')
		status, newmtype_excit = dict_sort_by_seq(mtype_excit, sq_ex[_i])
		if not status: exit('er2')
	

		mpl.rcParams['font.size'] = 33.0
		ax = plt.subplot(5, 3, (_i * 3) + 1, frameon=False)
		ax.barh(['_', layer, ''], [0, i/(i+e)  + e/(i+e), 0], height=0.4, color='blue')
		ax.barh(['_', layer, ''], [0, e/(i+e), 0], height=0.4, color='red')
		ax.set_yticklabels(['', layer + ' '*5, ''])
		ax.set_xticklabels([])
		mpl.rcParams['font.size'] = 23.0
		ax.text(0.85, 0, str(i), color='blue')
		ax.text(0.05, 0, str(e), color='red')
		ax.set_frame_on(False)
		ax.set_yticklabels(['', layer + ' '*7, ''])
		mpl.rcParams['font.size'] = 15.0
		ax = plt.subplot(5, 3, (_i * 3) + 2, frameon=False)
		ax.pie(list(newmtype_excit.values()), labels=list(newmtype_excit.keys()), startangle=angelex[_i], radius = 1.3, counterclock=False,  textprops={'size': 'smaller'})
		ax = plt.subplot(5, 3, (_i * 3) + 3, frameon=False)
		ax.pie(list(newmtype_inhib.values()), labels=list(newmtype_inhib.keys()), startangle=angelin[_i], counterclock=False, radius = 1.3,  textprops={'size': 'smaller'})
	plt.savefig('0.05Fig30.png')


def dict_sort_by_seq(d, s):
	newd = {}
	for _d in d.keys():
		if not _d in s:
			return False, 0
	for _s in s:
		newd[_s] = d[_s]
	return sum(newd.values()) == sum(d.values()), newd

def get_mtype(cellid):
	global should_mtypes
	for i in should_mtypes:
		if cellid.startswith(i):
			return True, i
	return False, -1
def get_reverse(d):
	l = list(d.keys())
	l.reverse()
	res = {}
	for i in l:
		res[i] = d[i]
	return res

def plot_synapses(res_folder='../CloudResultSynapse'):
	global should_mtypes, should_etypes

	inhib, excit = get_inh_exi()

	em = {i:0 for i in inhib.keys()}
	im = {i:0 for i in inhib.keys()}
	if True:
	# if False:
	# if not os.path.exists(opj(res_folder, 'Mpre_syn.json')):
	# if os.path.exists(opj(res_folder, 'Mpre_syn.json')):
		dict_connections = {}
		# file = 'synapses_0_map.json'	
		# file = 'synapses_1_map.json'
		file = 'synapses_map.json'
		# file = 'synapses_3_map.json'
		
		with open(opj(res_folder, file), 'r') as outfile:
			result = json.load(outfile)
		Mpre_syn={}
		Mpost_syn={}
		for k in excit.keys():
			Mpre_syn[k]=[]#how many uniq preneurons cell id in this mtype
			Mpost_syn[k]=[]

		r_result = {}
		all_m_pre = []
		lll = len(result.keys())
		for post_syn, pre_synaptics in result.items():
			print(lll)
			lll -= 1
			status, mtype_post_syn = get_mtype(post_syn)
			if not status:
				print(pre_synaptics, post_syn)
				
			uniq = list(np.unique(list(pre_synaptics.values())))

			# Mpre_syn[mtype_post_syn].extend(uniq)
			Mpre_syn[mtype_post_syn].append(len(uniq))

			dpre_syn = {pre_syn:0 for pre_syn in uniq}
			
			for num, pre_syn in pre_synaptics.items():
				dpre_syn[pre_syn] += 1

			for pre_syn in uniq:
				status, mtype_pre_syn = get_mtype(pre_syn)
				if not status:
					print(pre_synaptics, post_syn)
				
				if inhib[mtype_pre_syn] == 1:
					im[mtype_post_syn] += dpre_syn[pre_syn]
				else:
					em[mtype_post_syn] += dpre_syn[pre_syn]

				r_result.setdefault(pre_syn, [])
				r_result[pre_syn].append(post_syn)

		r_result[pre_syn].append(post_syn)
		ll = len(r_result)
		for pre_syn, list_post_syn in r_result.items():
			print(ll)
			ll -= 1
			nums_post = len(np.unique(list_post_syn))
			status, mtype_pre_syn = get_mtype(pre_syn)
			if not status:
				print(pre_synaptics, post_syn)

			Mpost_syn[mtype_pre_syn].append(nums_post)

		for k in excit.keys():
			Mpost_syn[k] = np.mean(Mpost_syn[k])
			Mpre_syn[k] = np.mean(Mpre_syn[k])

		with open(opj(res_folder, 'em.json'), 'w') as outfile:
			json.dump(em, outfile)
		with open(opj(res_folder, 'im.json'), 'w') as outfile:
			json.dump(im, outfile)

		with open(opj(res_folder, 'Mpre_syn.json'), 'w') as outfile:
			json.dump(Mpre_syn, outfile)
		with open(opj(res_folder, 'Mpost_syn.json'), 'w') as outfile:
			json.dump(Mpost_syn, outfile)
	else:
		# with open(opj(res_folder, 'em.json'), 'r') as outfile:
		# 	em = json.load(outfile)
		# with open(opj(res_folder, 'im.json'), 'r') as outfile:
		# 	im = json.load(outfile)

		# res_folder = '../Connectome'
		with open(opj(res_folder, 'Mpost_syn.json'), 'r') as outfile:
			Mpost_syn = json.load(outfile)
		with open(opj(res_folder, 'Mpre_syn.json'), 'r') as outfile:
			Mpre_syn = json.load(outfile)
	
	
	Mpre_syn = {k.replace('-','_'):v for k, v in Mpre_syn.items()}
	Mpost_syn = {k.replace('-','_'):v for k, v in Mpost_syn.items()}
	
	ipost, ipre = {}, {}
	epost, epre = {}, {}
	for k in excit.keys():
		if excit[k] != 0:
			epost[k], epre[k] = Mpost_syn[k], Mpre_syn[k]
		elif inhib[k] != 0:
			ipost[k], ipre[k] = Mpost_syn[k], Mpre_syn[k]

	epre = get_reverse(epre)
	epost = get_reverse(epost)
	ipre = get_reverse(ipre)
	ipost = get_reverse(ipost)

	for showup, showdown, name in [[True, False, '2up'], [True, True, '1all'],[False, True, '3down']]:

		
		mpl.rcParams['font.size'] = 23.0
		print(showup, showdown, name)
		if showup and showdown:
			plt.figure(figsize=(13, 22))
			up1, up2 = 221, 222
			d1,d2 = 223, 224
		elif showup:
			plt.figure(figsize=(13, 13))
			up1, up2 = 121, 122
		elif showdown:
			plt.figure(figsize=(13, 18))
			d1,d2 = 121, 122

		if showup:
			ax = plt.subplot(up1, frameon=False)
			ax.barh(list(epost.keys()), list(epre.values()), color='red')

			ax.set_xlim([800, 0])
			ax.set_yticklabels([])
			ax = plt.subplot(up2, frameon=False)
			ax.barh(list(epost.keys()), list(epost.values()), color='blue')
			ax.set_xlim([0, 800])

		if showdown:
			ax = plt.subplot(d1, frameon=False)
			ax.barh(list(ipost.keys()), list(ipre.values()), color='red')
			ax.set_xlim([500, 0])
			ax.set_yticklabels([])
			ax = plt.subplot(d2, frameon=False)
			ax.barh(list(ipost.keys()), list(ipost.values()), color='blue')
			ax.set_xlim([0, 500])
			print(sum(list(ipost.values()))+ sum(list(epost.values())), sum(list(ipre.values()))+ sum(list(epre.values())))
		
		# namefile = opj('11', name + '.png')
		namefile = opj('11', name + '_NEw2.png')
		plt.savefig(namefile)
		print(namefile)
		plt.cla()

	#Fig32
	plt.rcParams.update({'font.size': 10})
	for k in em.keys():
		im[k]=im[k]/(im[k] + em[k])
		em[k]=em[k]/(im[k] + em[k])
		em[k]+=im[k]
	mpl.rcParams['font.size'] = 15.0
	plt.figure(figsize=(20, 10))
	ax = plt.subplot(2,1,1, frameon=False)
	
	plt.bar(em.keys(), em.values(), color = 'red')
	plt.bar(im.keys(), im.values(), color = 'blue')
	# plt.xticks(range(len(im)), im.keys(), fontweight='heavy')
	plt.xticks(rotation=70)
	plt.ylim([0,1])
	# plt.show()
	plt.savefig('0.05fig32.png')

	


def pplot_per_layer(template_cell_ids, t, result, folder, figsize, t1=0,t2=0):

	template_cell_ids.reverse()
	if not (template_cell_ids[0].startswith('L6') and template_cell_ids[-1].startswith('L1')):
		print(template_cell_ids[0], template_cell_ids[-1])
		exit('-')

	mpl.rcParams['font.size'] = 15.0
	l = len(template_cell_ids)
	print('pplot_per_layer', len(template_cell_ids))
	print(template_cell_ids[0], template_cell_ids[-1])
	plt.figure(figsize=figsize)
	
	ix = []
	for i in range(len(template_cell_ids)):
		ix.append(i*100)

	for i, cell_id in enumerate(template_cell_ids):
		add = ix[i]
		newres = [i + add for i in result[cell_id]]
		plt.plot(t, newres, color='black')

	plt.xlabel('t')
	plt.yticks([i-50 for i in ix], template_cell_ids)
	plt.savefig(opj(folder, '1plot_per_layer.png'))


def pvlines_per_layer(template_cell_ids, t, result, folder, t1=0,t2=0, figsize=(13,9)):
	global should_mtypes
	
	template_cell_ids.reverse()
	if not (template_cell_ids[0].startswith('L6') and template_cell_ids[-1].startswith('L1')):
		print(template_cell_ids[0], template_cell_ids[-1])
		exit('-')
	mpl.rcParams['font.size'] = 18.0
	l = len(template_cell_ids)
	
	print('pvlines_per_layer', len(template_cell_ids))
	print(template_cell_ids[0], template_cell_ids[-1])
	cell_names = {'L1':0, 'L23':0, 'L4':0, 'L5':0, 'L6':0}

	for i in template_cell_ids:
		layer = i.split('_')[0]
		cell_names[layer] += 1
	label = []
	
	ix = []
	nms = 0
	for layer in ['L6','L5','L4','L23','L1']:
		v = int(cell_names[layer] / 2)
		nms += v
		ix.append(nms)
		label.append(layer)
		nms += v
	print(ix)
	print(label)

	plt.figure(figsize=figsize)
	for i, cell_id in enumerate(template_cell_ids):
		for spike in result[cell_id]:		
			plt.vlines(spike, i, i+1)

	plt.yticks(ix, label, size='small')
	plt.ylim([0,l+1])
	if t1 != 0 and t2 != 0:
		plt.xlim([t1,t2])
	# plt.xlim([0,2500])
	plt.xlabel('t')
	plt.savefig(opj(folder, '1pvline.png'))

def hihistogramm_layer(t1,t2,t, all_spike_histogramm, folder, figsize):
	if t1 !=0 and t2 != 0:
		b = np.linspace(t1,t2, 200, endpoint=True)
	else:
		b = np.linspace(0,max(t), 500, endpoint=True)


	spike_histogramm = {'L1':[], 'L23':[], 'L4':[], 'L5':[], 'L6':[]}
	for k,v in all_spike_histogramm.items():
		spike_histogramm[k.split('_')[0]].extend(v)
	l = len(spike_histogramm)
	mpl.rcParams['font.size'] = 15.0
	plt.figure(figsize=figsize)
	
	for i, (cell_id, m_spikes) in enumerate(spike_histogramm.items()):
		ax = plt.subplot(l, 1, i+1, frameon=False)
		plt.hist(m_spikes, b, color='black')

		plt.xlim([b[0],b[-1]])
		# plt.xlim([0,2500])
		plt.ylim([0,70])
		h = ax.set_ylabel(cell_id)
		h.set_rotation(0)
		ax.set_yticklabels([])
		if i < l - 1:
			ax.set_xticklabels([])

	h = ax.set_xlabel('t')
	plt.savefig(opj(folder, 'hihistogramm_layer.png'))


def hihistogramm_mtype(t1,t2,t, spike_histogramm, folder, figsize):
	if t1 !=0 and t2 != 0:
		b = np.linspace(t1,t2, 200, endpoint=True)
	else:
		b = np.linspace(0,max(t), 250, endpoint=True)

	l = len(spike_histogramm)
	plt.figure(figsize=figsize)
	mpl.rcParams['font.size'] = 15.0
	for i, (cell_id, m_spikes) in enumerate(spike_histogramm.items()):
		ax = plt.subplot(l, 1, i+1, frameon=False)
		if len(m_spikes):
			plt.hist(m_spikes, b, color='black')

		plt.xlim([b[0],b[-1]])
		# plt.xlim([0,2500])
		plt.ylim([0,50])
		
		h = ax.set_ylabel(cell_id)
		h.set_rotation(0)
		ax.set_yticklabels([])
		if i < l - 1:
			ax.set_xticklabels([])

	h = ax.set_xlabel('t')
	plt.savefig(opj(folder, 'hihistogramm_mtype.png'))

def plot_parallel_cores(t1_t2, folder='', res_folder='result'):
	global should_mtypes, should_etypes
	
	inhib, excit = get_inh_exi()
	create_fig30(inhib, excit, res_folder)
	exit()
	num_neurons = 0
	t = []
	short_soma_v = {}
	cell_spike_times = {}
	spike_times = {}
	hisst = {}
	
	for cell_id in should_mtypes:
		hisst[cell_id] = []
	layer_nums = {'L1':0, 'L23':0, 'L4':0, 'L5':0, 'L6':0}
	for file in os.listdir(res_folder):
		if file.startswith('host_') and file.endswith('.json'):
			print(file)

			with open(opj(res_folder, file), 'r') as outfile:
				result = json.load(outfile)
			if not len(t):
				t = result['t']

			for iii, (cell_id, soma_v) in enumerate(result['record_soma_v'].items()):
				spiket = result['record_spike_times'][cell_id]
				for m_cell_id in should_mtypes:
					if cell_id.startswith(m_cell_id):
						break
				
				layer_nums[cell_id.split('_')[0]] += 1
				
				spike_times.setdefault(m_cell_id, [])
				spike_times[m_cell_id].extend(spiket)

				cell_spike_times[cell_id] = np.asarray(spiket)

				isi = []
				if len(spiket) > 1:
					# isi = [r-l for r,l in zip(spiket[1:], spiket[:-1])]
					isi = spiket
					hisst[m_cell_id].extend(isi)
				
				soma_v = result['record_soma_v'][cell_id]


				soma_v = np.asarray(soma_v)
				
				if not m_cell_id in short_soma_v.keys():
					short_soma_v[m_cell_id] = soma_v
				# else:
				# 	was_soma_v = short_soma_v[m_cell_id]
				# 	if np.std(soma_v) > np.std(was_soma_v):
				# 		short_soma_v[m_cell_id] = soma_v
				
				num_neurons += 1
		# break
	for ix in range(len(t1_t2)):
		_t1, _t2 = t1_t2[ix]

		short_t = t
		frstt = -1
		lstt = -1
		if _t1 > 0 and _t2 > 0:
			for i in range(len(t)):
				if frstt == -1 and _t1 < t[i]:
					frstt = i
				if lstt == -1 and t[i] > _t2:
					lstt = i
					break
			short_t = [i for i in t if _t1 < i < _t2]

		ssshort_soma_v = short_soma_v.copy()
		if frstt != -1 and lstt != -1:
			for k, v in short_soma_v.items():
				ssshort_soma_v[k] = short_soma_v[k][frstt:lstt]
		
		sshort_spike_times = spike_times.copy()
		if _t1 > 0 and _t2 > 0:
			for k in sshort_spike_times.keys():
				sshort_spike_times[k] = [i for i in sshort_spike_times[k] if _t1 < i < _t2]

		short_cell_spike_times = cell_spike_times.copy()
		if _t1 > 0 and _t2 > 0:
			for k in short_cell_spike_times.keys():
				short_cell_spike_times[k] = np.asarray([i for i in short_cell_spike_times[k] if _t1 < i < _t2])
		short_hisst = hisst.copy()
		if _t1 > 0 and _t2 > 0:
			for k in short_hisst.keys():
				short_hisst[k] = [i for i in short_hisst[k] if _t1 < i < _t2]
		short_folder = ''
		if folder != '':
			short_folder = folder
			if _t1 != 0 or _t2 != 0:
				short_folder = f'{folder}_{_t1}-{_t2}'
			shutil.rmtree(short_folder, ignore_errors=True)
			os.makedirs(short_folder, exist_ok=True)
		
		plot_result(short_t, _t1,_t2, short_cell_spike_times, layer_nums, ssshort_soma_v, short_hisst, sshort_spike_times, short_folder, res_folder, num_neurons)
	
def plot_result(t, t1,t2, cell_spike_times, layer_nums, short_soma_v, hisst, spike_times, folder, res_folder, num_neurons):
	print('\n\n\n\n')
	print(folder)
	print(t1, t2, 't1, t2')
	print(layer_nums)
	print(len(cell_spike_times.keys()), num_neurons, '----------------')

	pattern_dict = {}
	for cell_id in spike_times.keys():
		spike_times[cell_id] = np.asarray(spike_times[cell_id])
	step = 10
	if t1 != 0 and t2 != 0:
		for window in range(t1-10, t2, step):
			layer_activity = {'L1':0, 'L23':0, 'L4':0, 'L5':0, 'L6':0}
			for cell_id, spiket in spike_times.items():
				layer = cell_id.split('_')[0]
				layer_activity[layer] += len(spiket[(spiket > window) & (spiket < window + 10)])
			if sum(layer_activity.values())>0:
				layer_activity = dict(sorted(layer_activity.items(), key=lambda x: x[1], reverse=True))
				# print(f'{window}:{window + step} 	{layer_activity}')
				print(f'{window}:{layer_activity}')
				pattern_dict[window] = layer_activity
		print(pattern_dict)
		bursting_per_time(layer_nums, pattern_dict, folder)
	sorted_cell_spike, wascellid = get_sorted_cellids(list(cell_spike_times.keys()), should_mtypes)
	sorted_short_soma, wassoma = get_sorted_cellids(list(short_soma_v.keys()), should_mtypes)
	if not os.path.exists(opj(res_folder, 'cell.json')):
		result = {'cell_ids':sorted_cell_spike}
		with open(opj(res_folder, 'cell.json'), 'w') as fd:
			json.dump(result, fd)

	print('nums:', wascellid, wascellid==len(sorted_cell_spike)==len(np.unique(sorted_cell_spike)),  'num_neurons:', num_neurons)
	if folder == '':
		exit('3-3--3-3-')
	figsize=(7,14)
	if t1 == 0 and t2 == 0:
		figsize=(14,14)

	hihistogramm_mtype(t1,t2,t, hisst, folder=folder, figsize=figsize)
	hihistogramm_layer(t1,t2,t, hisst, folder=folder, figsize=figsize)
	
	pvlines_per_layer(sorted_cell_spike, t, cell_spike_times, folder=folder, t1=t1,t2=t2, figsize=figsize)
	pplot_per_layer(sorted_short_soma, t, short_soma_v, folder=folder, t1=t1,t2=t2, figsize=figsize)

	


if __name__ == '__main__':
	opj = os.path.join
	should_mtypes = ['L1_DAC','L1_NGC_DA','L1_NGC_SA','L1_HAC','L1_DLAC','L1_SLAC','L23_PC','L23_MC','L23_BTC','L23_DBC','L23_BP','L23_NGC','L23_LBC','L23_NBC','L23_SBC','L23_ChC','L4_PC','L4_SP','L4_SS','L4_MC','L4_BTC','L4_DBC','L4_BP','L4_NGC','L4_LBC','L4_NBC','L4_SBC','L4_ChC','L5_TTPC1','L5_TTPC2','L5_UTPC','L5_STPC','L5_MC','L5_BTC','L5_DBC','L5_BP','L5_NGC','L5_LBC','L5_NBC','L5_SBC','L5_ChC','L6_TPC_L1','L6_TPC_L4','L6_UTPC','L6_IPC','L6_BPC','L6_MC','L6_BTC','L6_DBC','L6_BP','L6_NGC','L6_LBC','L6_NBC','L6_SBC','L6_ChC']
	should_etypes = ['cAC', 'bAC', 'cNAC', 'bNAC', 'dNAC', 'cSTUT', 'bSTUT', 'dSTUT', 'cIR' , 'bIR', 'cADpyr']
	

	plot_synapses(res_folder='result/')

	t1_t2 = [[0,0], [210, 310], [1450, 1550], [2870, 2970]]#1442
	plot_parallel_cores(t1_t2, folder='plots', res_folder='result/')
