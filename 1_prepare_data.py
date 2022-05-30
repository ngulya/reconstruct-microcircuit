import os
import matplotlib
import numpy as np
import sys
import pylab
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import json
import shutil
import fileinput
from filecmp import dircmp
import copy
import time
from scipy.spatial import distance
import random as rnd
import re
rnd.seed(13)


def get_neurons_path(circuit_folder):
	l_neuron_folders = []
	for neuron_folders in os.listdir(circuit_folder):
		if os.path.isdir(opj(circuit_folder, neuron_folders)):
			l_neuron_folders.append(neuron_folders)
	return l_neuron_folders


def create_Mpost_pre_syn_json_2022(anatomy_path):
	Mpre_syn_nums, Mpost_syn_nums = {}, {}
	with open(anatomy_path, 'r') as fd:
		anatomy = json.load(fd)

	
	for connection, v in anatomy.items():
		pre, post = connection.split(':')
		pre, post = pre.replace('-', '_'), post.replace('-', '_')
		Mpre_syn_nums.setdefault(pre, 0)
		Mpost_syn_nums.setdefault(post, 0)

		Mpre_syn_nums[pre] += v['number_of_divergent_neuron_mean']
		Mpost_syn_nums[post] += v['number_of_convergent_neuron_mean']
	return Mpre_syn_nums, Mpost_syn_nums

def create_synapse_map_v2(anatomy_path, folder=''):
	
	# with open(opj(folder,'Mpost_syn.json'), 'r') as outfile:
	# 	Mpost_syn_nums = json.load(outfile)
	# with open(opj(folder,'Mpre_syn.json'), 'r') as outfile:
	# 	Mpre_syn_nums = json.load(outfile)
	# Mpre_syn_nums = {k.replace('-','_'):v for k, v in Mpre_syn_nums.items()}
	# Mpost_syn_nums = {k.replace('-','_'):v for k, v in Mpost_syn_nums.items()}
	#?????
	Mpre_syn_nums, Mpost_syn_nums = create_Mpost_pre_syn_json_2022(anatomy_path)
	# print(Mpost_syn_nums['L1_SLAC'])
	# print(Mpost_syn_nums['L1_DLAC'])
	# print(Mpost_syn_nums['L1_DAC'])
	# print(Mpost_syn_nums['L1_SAC'])
	
	# exit()
	#difference num pre post
	inhib, excit = get_inh_exi()

	str_id_n_id, n_id_str_id, gid_syn_pre_cell_type, gid_R_syn_pre_cell_type, gid_num_pre_syn_mtype = get_lists()
	


	
	with open(anatomy_path, 'r') as fd:
		anatomy = json.load(fd)

	synapses_map = {}
	used_Mpost_syn_nums, used_Mpre_syn_nums = {},{}
	for cell_id_str in n_id_str_id:
		m_type, _ = get_mtype(cell_id_str, inhib)
		used_Mpost_syn_nums.setdefault(m_type, {})
		used_Mpre_syn_nums.setdefault(m_type, {})
		used_Mpost_syn_nums[m_type][cell_id_str] = 0
		used_Mpre_syn_nums[m_type][cell_id_str] = 0
	
	
	## str_id_n_id[post_nrn_str_id] = len(n_id_str_id) # {'L4_BP_bAC217_2': 0}
	## n_id_str_id[index] = post_nrn_str_id ['L4_BP_bAC217_2']
	## gid_syn_pre_cell_type = [0]{0:'L1_DAC', 1:'L4_SP'} [index from all 30000*multiplier neurons][synapse_id:mtype_pre]
	## gid_R_syn_pre_cell_type = [0]{'L1_DAC':0, 'L4_SP':1} [index from all 30000*multiplier neurons][mtype_pre:synapse_id]
	## gid_num_pre_syn_mtype = [0]{'L1_DAC':31, 'L4_SP':20} [index from all 30000*multiplier neurons][mtype_pre:num] | number of pre neuron mtypes

	'''
	used_Mpre_syn_nums
	{mtype:{mtype_etype_dplc_IX:0}}
	{....'L6_BTC': {'L6_BTC_cNAC_1_dplc_0': 0, 'L6_BTC_bAC_1_dplc_0': 0, 'L6_BTC_cAC_1_dplc_0': 0}.....}
	'''
	
	for mtype in should_mtypes:

		should_be_post = max([1, int(Mpost_syn_nums[m_type]/len(used_Mpost_syn_nums[m_type].keys()))])
		should_be_pre = max([1,int(Mpre_syn_nums[m_type]/len(used_Mpre_syn_nums[m_type].keys()))])
		for k in used_Mpost_syn_nums[m_type].keys():
			used_Mpost_syn_nums[m_type][k] = should_be_post
		for k in used_Mpre_syn_nums[m_type].keys():
			used_Mpre_syn_nums[m_type][k] = should_be_pre


	# print(excit['L1_SAC'])
	# print(excit['L1_DAC'])
	# exit()

	#used_Mpost_syn_nums -> {'L1_DAC':{'L1_DAC_cUdf_13':100, 'L1_DAC_cUdf_113':100}}
	ll = len(n_id_str_id)
	for _i, post_syn_num in enumerate(range(ll)):
		post_syn_str_id = n_id_str_id[post_syn_num]
		print('create_synapse_map_v2: ', ll-_i, post_syn_str_id)
		post_syn_mtype, this_inhib = get_mtype(post_syn_str_id, inhib)

		#num of pre mtypes
		pre_syn_mtype_num = gid_num_pre_syn_mtype[post_syn_num].copy()#{'L4_PC': 7, 'L4_SP': 2, 'L4_BTC': 3, 'L6_LBC': 2, 'L23_DBC': 3, 'L5_MC': 8, 'L5_NBC': 6}
		R_pre_syn_mtype = gid_R_syn_pre_cell_type[post_syn_num].copy()#{'L4_PC': [0, 1, 6, 10, 11, 17, 25], 'L4_SP': [2, 3]}....

		synapses_map[post_syn_str_id] = gid_syn_pre_cell_type[post_syn_num].copy()#{0: 'L4_PC', 1: 'L4_PC', 2: 'L4_SP', 3: 'L4_SP', 4: 'L4_SP', 5: 'L4_SP', .....

		synaptic_connections = sum(pre_syn_mtype_num.values())
		pre_syn_uniq_cell = Mpre_syn_nums[post_syn_mtype]

		map_pre_syn_uniq_cell = pre_syn_mtype_num.copy()
		if pre_syn_uniq_cell < synaptic_connections:
			for pre_syn_mtype, num_synapses in pre_syn_mtype_num.items():####z1
				_pre_syn_mtype = pre_syn_mtype
				if pre_syn_mtype in ['L1_NGC_DA', 'L1_NGC_SA']:
					_pre_syn_mtype = pre_syn_mtype.split('_')
					_pre_syn_mtype = _pre_syn_mtype[0] + '_' + _pre_syn_mtype[1] + '-' + _pre_syn_mtype[2]
				_post_syn_mtype = post_syn_mtype
				if post_syn_mtype in ['L1_NGC_DA', 'L1_NGC_SA']:
					_post_syn_mtype = post_syn_mtype.split('_')
					_post_syn_mtype = _post_syn_mtype[0] + '_' + _post_syn_mtype[1] + '-' + _post_syn_mtype[2]
				name = _pre_syn_mtype+':'+_post_syn_mtype
				if name in anatomy:
					mean_number_of_synapse_per_connection = anatomy[name]['mean_number_of_synapse_per_connection']
				else:
					mean_number_of_synapse_per_connection = 1
				map_pre_syn_uniq_cell[pre_syn_mtype] = max([1, int(num_synapses/mean_number_of_synapse_per_connection)])

			kof = pre_syn_uniq_cell/sum(map_pre_syn_uniq_cell.values())####z2
			map_pre_syn_uniq_cell = {k:max([1, int(v*kof)]) for k, v in map_pre_syn_uniq_cell.items()}
			##{'L4_PC': 1, 'L4_SP': 1, 'L4_BTC': 3, 'L6_LBC': 1, 'L23_DBC': 1, 'L5_MC': 2, 'L5_NBC': 1}#cells
		for pre_syn_mtype, num_synapses in pre_syn_mtype_num.items():#{'L4_PC': 7, 'L4_SP': 2, 'L4_BTC': 3, 'L6_LBC': 2, 'L23_DBC': 3, 'L5_MC': 8, 'L5_NBC': 6}
			# available_gid = has_this_mtype(n_id_str_id, post_syn_num, pre_syn_mtype, map_pre_syn_uniq_cell[pre_syn_mtype], used_Mpost_syn_nums[pre_syn_mtype].copy())
			# print(post_syn_str_id, pre_syn_mtype)

			available_gid = has_this_mtype_z3(post_syn_str_id, str_id_n_id, map_pre_syn_uniq_cell[pre_syn_mtype], used_Mpost_syn_nums[pre_syn_mtype].copy())
			num_available = len(available_gid)

			if num_available:
				drop = False
				if num_synapses <= num_available:
					drop = True
				used_cell_ids = {}
				for syn_id in R_pre_syn_mtype[pre_syn_mtype]:
					# print(drop, post_syn_mtype)
					cell_id = rnd.choice(available_gid)
					if drop:
						available_gid.remove(cell_id)
					used_cell_ids[n_id_str_id[cell_id]] = 0
					synapses_map[post_syn_str_id][syn_id] = n_id_str_id[cell_id]
				for used in used_cell_ids:
					used_Mpost_syn_nums[pre_syn_mtype][used] -= 1
			else:
				for syn_id in R_pre_syn_mtype[pre_syn_mtype]:
					synapses_map[post_syn_str_id][syn_id] = 'no-data'
					exit('Error 1')
	if not os.path.exists('result'):
		os.makedirs('result')
	with open(f'result/synapses_map.json', 'w') as outfile:
		json.dump(synapses_map, outfile)
	with open(f'result/used_Mpost_syn_nums.json', 'w') as outfile:
		json.dump(used_Mpost_syn_nums, outfile)
	return synapses_map


def change_str_in_hoc(filename, text_to_search, text_to_replace):
	with fileinput.FileInput(filename, inplace=True) as file:
		for line in file:
			print(line.replace(text_to_search, text_to_replace), end='')
	return 1

def change_template_name(filename, old, new, has_old_name=True):
	core_word = 'begintemplate'
	old_was = ''
	with fileinput.FileInput(filename, inplace=True) as file:
		for line in file:
			if line.find(core_word) >=0:
				if has_old_name and line.find(old) < 0:
					exit(f'change_template_name({filename}, {old}, {new})')
				old_was	= line.split()[1]
				line = f'{core_word} {new}\n'
				core_word = 'endtemplate'
			print(line, end='')
	return old_was

def add_new_objects(filename, new_obj):
	core_word = 'begintemplate'
	changed = ''
	with fileinput.FileInput(filename, inplace=True) as file:
		for line in file:
			if line.find(core_word) >=0:
				line = f'{line}  public {new_obj}\n  objref {new_obj}\n'
				changed = line
			print(line, end='')
	return changed

def keep_axon(filename):
	first_replace = True
	with fileinput.FileInput(filename, inplace=True) as file:
		for line in file:
			if first_replace and line.find('replace_axon()') != -1:
				line = line.replace('replace_axon', '//replace_axon')
				first_replace = False
			print(line, end='')
	return 1

def load_neurons_data(d_map=None, axon_use=True, circuit_folders=['MyLayers'], name_cell=None):
	
	file_loc = os.path.dirname(os.path.abspath(__file__))
	
	for folder in ['synapses', 'morphology', 'biophysics', 'template', 'current_amps']:
		shutil.rmtree(folder, ignore_errors=True)
		os.makedirs(folder, exist_ok=True)
	shutil.rmtree('mechanisms', ignore_errors=True)
	
	have_mechanisms = False
	neuron_parent_folders = []
	for folder in circuit_folders:
		neuron_parent_folders.append(get_neurons_path(folder))

	for _num_fold, neuron_folders in enumerate(neuron_parent_folders):
		for cell_id in neuron_folders:
			if not name_cell is None and cell_id != name_cell:
				continue
			print(cell_id)
			folder = opj(circuit_folders[_num_fold], cell_id)
			with open(opj(folder, 'cellinfo.json'), 'r') as fd:
				cell_info = json.load(fd)
			cell_id = cell_id.replace('-', '_')
			biophysic_old_name = cell_info["cell name"].split('_')[0] + '_biophys'
			biophysic_new_name = cell_id + '_biophys'
			synapse_new_name = 'synapses_' + cell_id
			morphology_new_name = 'morphology_' + cell_id

			biophysics_hoc_old = opj(folder, 'biophysics.hoc')
			biophysics_hoc_new = opj('biophysics', f'{cell_id}.hoc')
			shutil.copy(biophysics_hoc_old, biophysics_hoc_new)
			change_template_name(biophysics_hoc_new, biophysic_old_name, biophysic_new_name)


			current_amps_old = opj(folder, 'current_amps.dat')
			current_amps_new = opj('current_amps', f'{cell_id}.dat')
			shutil.copy(current_amps_old, current_amps_new)

			morphology_old = opj(folder, 'morphology')
			morphology_new = opj('morphology', cell_id)
			shutil.copytree(morphology_old, morphology_new)

			morphology_hoc_old = opj(folder, 'morphology.hoc')
			morphology_hoc_new = opj(morphology_new, 'morphology.hoc')
			shutil.copy(morphology_hoc_old, morphology_hoc_new)
			change_str_in_hoc(opj(morphology_new, 'morphology.hoc'), 'morphology/', opj(file_loc, morphology_new+'/'))
			morphology_old_name = change_template_name(morphology_hoc_new, None, morphology_new_name, False)

			synapses_old = opj(folder, 'synapses')
			synapses_new = opj('synapses', cell_id)
			shutil.copytree(synapses_old, synapses_new)
			
			change_str_in_hoc(opj(synapses_new, 'synapses.hoc'), 'synapses/', opj(file_loc, synapses_new+'/'))
			change_str_in_hoc(opj(synapses_new, 'synapses.hoc'), 'synapse_list.append(synapse)', 'synapse.synapseID = synapse_id\nsynapse_list.append(synapse)')
			synapse_old_name = change_template_name(opj(synapses_new, 'synapses.hoc'), None, synapse_new_name, False)
			
			template_hoc_old = opj(folder, 'template.hoc')
			template_hoc_new = opj('template', f'{cell_id}.hoc')

			shutil.copy(template_hoc_old, template_hoc_new)
			change_str_in_hoc(template_hoc_new, 'morphology.hoc', opj(file_loc, morphology_hoc_new))
			change_str_in_hoc(template_hoc_new, 'biophysics.hoc', opj(file_loc, biophysics_hoc_new))
			change_str_in_hoc(template_hoc_new, 'synapses/synapses.hoc', opj(file_loc, synapses_new, 'synapses.hoc'))

			#change template names in addition mechs
			change_str_in_hoc(template_hoc_new, biophysic_old_name, biophysic_new_name)
			change_str_in_hoc(template_hoc_new, synapse_old_name, synapse_new_name)
			change_str_in_hoc(template_hoc_new, morphology_old_name, morphology_new_name)
			change_template_name(template_hoc_new, None, cell_id, False)
			# add_new_objects(template_hoc_new, 'soma_v')
			if axon_use:
				keep_axon(template_hoc_new)
			if not have_mechanisms:
				shutil.copytree(opj(folder, 'mechanisms'), 'mechanisms')
				have_mechanisms = True
			else:
				new_files = dircmp(opj(folder, 'mechanisms'), 'mechanisms').left_only
				for file in new_files:
					print('\n\n\n\n\n\n\nnrnivmodl mechanisms')
					shutil.copy(opj(folder, 'mechanisms', file), opj('mechanisms', file))

	return 1

def init_weights(mtype, etype, cell_mtype, cell_etype):
	EMtype, MEtype = {}, {}
	numx = 0
	for i in range(len(cell_mtype)):
		e = cell_etype[i]
		m = cell_mtype[i]
		MEtype.setdefault(m, {})
		EMtype.setdefault(e, {})

		EMtype[e][m] = mtype[m]
		MEtype[m][e] = etype[e]
		numx += 1
	return MEtype, EMtype


def get_equal_sum_mtype_etype(etype, mtype):
	nume, numm = sum(etype.values()), sum(mtype.values())

	while nume - numm != 0:
		e = rnd.choice(list(etype.keys()))
		if nume > numm:
			if etype[e] > 1:
				etype[e] -= 1
		else:
			etype[e] += 1
		nume = sum(etype.values())
	return etype, mtype

def change_layer_cell_number(lst_layer_cell, multiplier, etype, mtype):

	multiplier_me = sum(mtype.values())/sum(etype.values())
	for _type in etype.keys():
		etype[_type] = int(multiplier_me * etype[_type])
	

	imtype, ietype = mtype.copy(), etype.copy()

	for m in mtype:
		mtype[m] = max([int(mtype[m]*multiplier), 1])
	for e in etype:
		etype[e] = max([int(etype[e]*multiplier), 1])


	cells = []

	cell_mtype, cell_etype = [], []
	
	for cell in lst_layer_cell:#mtype only
		cw = cell
		cell = cell.split('_')
		_etype = cell[-1]
		cell = '_'.join(cell[:-1])

		if _etype in etype and cell in mtype:
			cell_mtype.append(cell)
			cell_etype.append(_etype)
			cells.append(cw)

	# print('\ncells:', cells)	
	#### for metype in cells:#['L1_NGC-SA', 'L1_LAC', 'L1_DAC', 'L1_NGC-DA', 'L1_HAC', 'L1_SLAC', 'L1_DLAC', 'L1_SAC']
	#### 	print('metype', metype)
	#### 	for _mtype in mtype.keys():#{'L1_DAC': 5, 'L1_DLAC': 2, 'L1_HAC': 9, 'L1_NGC-SA': 5, 'L1_NGC-DA': 7, 'L1_SLAC': 4}			
	#### 		if _mtype in metype:
	#### 			_etype_with_shit = metype.replace(_mtype+'_', '')
	#### 			print('	', _mtype,_etype_with_shit)
	#### 			if _mtype != _etype_with_shit:
	#### 				for _etype in etype.keys():#{'cAC': 1, 'cSTUT': 1, 'cNAC': 24, 'bNAC': 6, 'cIR': 1}
	#### 					if _etype in _etype_with_shit:
	#### 						cell_etype.append(_etype)
	#### 						cell_mtype.append(_mtype)
	#### 						break

	etype, mtype = get_equal_sum_mtype_etype(etype, mtype)


	MEtype, EMtype = init_weights(mtype, etype, cell_mtype, cell_etype)
	for m, e in MEtype.items():
		num_m = mtype[m]
		sum_e_types = sum(e.values())
		# print(m, e, 'sum_e_types',sum_e_types, num_m,'num_m')
		for _e, _num in e.items():
			v = max((int(_num*num_m/sum_e_types), 1))
			e[_e], EMtype[_e][m] = v, v
		# print(mtype[m],'	|',m, '	', e)#list(e.keys()), sum(e.values()))

	# for e, m in EMtype.items():
	# 	print(etype[e],'	|',e, '	', m)#ist(m.keys()), sum(m.values()))

	newmtype, newetype = get_new_statistics(mtype, etype, MEtype)
	# print('---')
	# print(imtype, 'imtype real')
	# print(mtype, 'imtype real normalized')
	# print(newmtype, 'imtype new')
	# print('---')
	# print(ietype, 'ietype real')
	# print(etype, 'ietype real normalized')
	# print(newetype, 'ietype new')
	
	return get_new_duplicated_list(MEtype, cells, cell_mtype, cell_etype)	

def get_cell(m, e, cells, cell_mtype, cell_etype):
	for i in range(len(cells)):# ['L1_NGC-SA', 'L1_LAC', 'L1_DAC', 'L1_NGC-DA', 'L1_HAC', 'L1_SLAC', 'L1_DLAC', 'L1_SAC']

		if cell_mtype[i] == m and cell_etype[i] == e:
			#['cNAC', 'cNAC', 'cAC', 'cNAC', 'bNAC', 'bNAC', 'cAC', 'cIR', 'cSTUT', 'bNAC', 'cNAC', 'cNAC', 'cNAC', 'bNAC']
			# ['L1_NGC-SA', 'L1_DAC', 'L1_NGC-DA', 'L1_HAC', 'L1_DAC', 'L1_HAC', 'L1_SLAC', 'L1_HAC', 'L1_NGC-DA', 'L1_SLAC', 'L1_SLAC', 'L1_DLAC', 'L1_NGC-DA', 'L1_NGC-DA']
			return cells[i]

def get_new_duplicated_list(MEtype, cells, cell_mtype, cell_etype):
	resulting = {}
	for m, _e in MEtype.items():
		# print(m, _e)
		for e, nums in _e.items():
			cell_id = get_cell(m, e, cells, cell_mtype, cell_etype)
			# print('	', e, nums, '-> ',cell_id, ':', nums)

			if cell_id in resulting.keys():
				exit('cell_id was in resulting')
			resulting[cell_id] = nums
	return resulting
	
def get_new_statistics(mtype, etype, MEtype):
	newmtype, newetype = mtype.copy(), etype.copy()
	for k in mtype.keys():
		newmtype[k] = 0

	for k in etype.keys():
		newetype[k] = 0

	for m, e in MEtype.items():
		newmtype[m] = sum(e.values())
		for _e, _n in e.items():
			newetype[_e] += _n
	return newmtype, newetype

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
		if inhib[k] != 0 and excit[k] != 0:
			print(k,inhib[k], excit[k])
		if inhib[k] > excit[k]:
			inhib[k] = 1
			excit[k] = 0
		else:
			inhib[k] = 0
			excit[k] = 1
	return inhib, excit

def duplicated_map(multiplier, json_folder='BBPjson', folder_layers='AllLayers'):
	'''
	multiplier = 0.1 means that only 10% from all 30k neurons will be modeleted
	'''
	# print('duplicated_map')
	cells = {'L1':[], 'L23':[],'L4':[],'L5':[],'L6':[]}
	for layer in os.listdir(folder_layers):
		layer_cell = opj(folder_layers, layer)
		if os.path.isdir(layer_cell):
			for cell in os.listdir(layer_cell):
				if os.path.isdir(opj(layer_cell, cell)):
					cells[layer].append(cell)

	# exit()
	# print('cells')
	# print(cells)

	with open(opj(json_folder, 'circuit_download.json'), 'r') as f:
		circuit_download = json.load(f)
	
	# print('circuit_download')
	# print(circuit_download)
	with open(opj(json_folder, 'layer_download.json'), 'r') as f:
		layer_download = json.load(f)
	
	num_cells, num_cells_in_crc = 0, 0
	# print(4873 +6081 +6453 +10080 + 144)
	d_map = {}
	for layer, lst_layer_cell in cells.items():

		_num_cells_in_crc = circuit_download["No. of neurons per layer"][layer]
		

		etype = layer_download[layer]["No. of neurons per electrical types"]
		mtype = layer_download[layer]["No. of neurons per morphological types"]
		# print('_num_cells_in_crc:', _num_cells_in_crc, sum(mtype.values()), 'mtype')
		# print('------------------------')
		d_map[layer] = change_layer_cell_number(lst_layer_cell, multiplier, etype, mtype)
		# print('------------------------')
		# print(d_map[layer])
		# exit()
		num_cells += sum(d_map[layer].values())
	print('will  be', num_cells)
	# exit('--')
	return d_map

def load_duplicated_neurons_data(d_map, axon_use=False, circuit_folders='AllLayers'):

	# print('load_duplicated_neurons_data')
	file_loc = os.path.dirname(os.path.abspath(__file__))
	
	for folder in ['synapses', 'morphology', 'biophysics', 'template', 'current_amps']:
		shutil.rmtree(folder, ignore_errors=True)
		os.makedirs(folder, exist_ok=True)
	shutil.rmtree('mechanisms', ignore_errors=True)
	

	'''
	#L1_SLAC_* no info in folder
	#L1_SAC_* no in anatomy, and zip in SAC startswith SLAC, so copy all to SLAC_*
	#L1_DLAC_cNAC no info in folder
	#L1_DAC_cNAC have info so copy zip with new name DLAC and replace in file all DLAC to DAC!
	'''
	for cell_folders in os.listdir(opj(circuit_folders, 'L1')):
		if cell_folders.startswith('L1_SAC') or cell_folders.startswith('L1_DAC_cNAC'):
			src_path = opj(circuit_folders, 'L1', cell_folders)
			dst_path = opj(circuit_folders, 'L1', cell_folders.replace('SAC', 'SLAC') if cell_folders.startswith('L1_SAC') else cell_folders.replace('DAC', 'DLAC'))
			slac_etype_files = os.listdir(dst_path)
			for files in os.listdir(src_path):
				if not files in slac_etype_files:
					# print(files, src_path)
					old = opj(src_path, files)
					new = opj(dst_path, files)
					if not os.path.isdir(old):
						shutil.copy(old, new)
						print(cell_folders, '|', files, '|', old, '->', 'new')

	have_mechanisms = False
	postfixs = ['_1','_2','_3','_4','_5']
	have_problem_with_cell = []
	numss = 0
	for layer in d_map.keys():
		numss += sum(d_map[layer].values())
	for layer in d_map.keys():
		for basic_cell_id, nums in d_map[layer].items():
			print('load_duplicated_neurons_data:', numss, basic_cell_id, nums)
			numss -= nums
			basic_cell_ids = []
			for _num_cpy in range(nums):

				if _num_cpy < 5:
					pfx = postfixs[_num_cpy]
				else:
					pfx = rnd.choice(postfixs)
					
				cell_folder = opj(circuit_folders, layer, basic_cell_id)
				folder = opj(cell_folder, basic_cell_id + pfx)

				if not os.path.exists(folder):#we should unzip L6_BPC_cADpyr231_2.zip -> L6_BPC_cADpyr_2 folder
					for filename in os.listdir(cell_folder):
						if filename.find('.zip') > 0:
							_filename = re.search(f'{basic_cell_id}.*.zip', filename)
							DLAC = False
							if basic_cell_id == 'L1_DLAC_cNAC':
								_filename = re.search(f'L1_DAC_cNAC.*.zip', filename)
								DLAC = True
							if not _filename is None:
								_filename = _filename.group()
								full_zip_path = opj(cell_folder, _filename)
								_filename = _filename.replace('.zip', '')
								will_be = basic_cell_id+'_'+_filename[-1]
								shutil.unpack_archive(full_zip_path, cell_folder)
								shutil.move(opj(cell_folder, _filename), opj(cell_folder, will_be))
								if DLAC:
									for file in ['cellinfo.json', 'createsimulation.hoc', 'run_RmpRiTau.py', 'template.hoc', 'run.py']:
										print(opj(cell_folder, will_be, file), 'DAC->DLAC')
										change_str_in_hoc(opj(cell_folder, will_be, file), 'DAC', 'DLAC')

				if not os.path.exists(opj(folder, 'cellinfo.json')):
					have_problem_with_cell.append(opj(folder, 'cellinfo.json'))
					break
				with open(opj(folder, 'cellinfo.json'), 'r') as fd:
					cell_info = json.load(fd)

				cell_id = basic_cell_id.replace('-', '_') + pfx
				cell_id = f'{cell_id}_dplc_{_num_cpy}'
				biophysic_old_name = cell_info["cell name"].split('_')[0] + '_biophys'
				biophysic_new_name = cell_id + '_biophys'
				synapse_new_name = 'synapses_' + cell_id
				morphology_new_name = 'morphology_' + cell_id

				biophysics_hoc_old = opj(folder, 'biophysics.hoc')
				biophysics_hoc_new = opj('biophysics', f'{cell_id}.hoc')
				shutil.copy(biophysics_hoc_old, biophysics_hoc_new)
				change_template_name(biophysics_hoc_new, biophysic_old_name, biophysic_new_name)

				current_amps_old = opj(folder, 'current_amps.dat')
				current_amps_new = opj('current_amps', f'{cell_id}.dat')
				shutil.copy(current_amps_old, current_amps_new)


				morphology_old = opj(folder, 'morphology')
				morphology_new = opj('morphology', cell_id)
				shutil.copytree(morphology_old, morphology_new)

				morphology_hoc_old = opj(folder, 'morphology.hoc')
				morphology_hoc_new = opj(morphology_new, 'morphology.hoc')
				shutil.copy(morphology_hoc_old, morphology_hoc_new)
				change_str_in_hoc(opj(morphology_new, 'morphology.hoc'), 'morphology/', opj(file_loc, morphology_new+'/'))
				morphology_old_name = change_template_name(morphology_hoc_new, None, morphology_new_name, False)

				synapses_old = opj(folder, 'synapses')
				synapses_new = opj('synapses', cell_id)
				shutil.copytree(synapses_old, synapses_new)
				
				change_str_in_hoc(opj(synapses_new, 'synapses.hoc'), 'synapses/', opj(file_loc, synapses_new+'/'))
				change_str_in_hoc(opj(synapses_new, 'synapses.hoc'), 'synapse_list.append(synapse)', 'synapse.synapseID = synapse_id\nsynapse_list.append(synapse)')
				synapse_old_name = change_template_name(opj(synapses_new, 'synapses.hoc'), None, synapse_new_name, False)
				
				template_hoc_old = opj(folder, 'template.hoc')
				template_hoc_new = opj('template', f'{cell_id}.hoc')
				shutil.copy(template_hoc_old, template_hoc_new)
				change_str_in_hoc(template_hoc_new, 'morphology.hoc', opj(file_loc, morphology_hoc_new))
				change_str_in_hoc(template_hoc_new, 'biophysics.hoc', opj(file_loc, biophysics_hoc_new))
				change_str_in_hoc(template_hoc_new, 'synapses/synapses.hoc', opj(file_loc, synapses_new, 'synapses.hoc'))

				#change template names in addition mechs
				change_str_in_hoc(template_hoc_new, biophysic_old_name, biophysic_new_name)
				change_str_in_hoc(template_hoc_new, synapse_old_name, synapse_new_name)
				change_str_in_hoc(template_hoc_new, morphology_old_name, morphology_new_name)
				change_template_name(template_hoc_new, None, cell_id, False)
				# add_new_objects(template_hoc_new, 'soma_v')
				
				if axon_use:
					keep_axon(template_hoc_new)
				
				if not have_mechanisms:
					shutil.copytree(opj(folder, 'mechanisms'), 'mechanisms')
					have_mechanisms = True
				else:
					new_files = dircmp(opj(folder, 'mechanisms'), 'mechanisms').left_only
					for file in new_files:
						shutil.copy(opj(folder, 'mechanisms', file), opj('mechanisms', file))


	if len(have_problem_with_cell) != 0:
		print('FUCK, Have a problem with some cells')
		for i in have_problem_with_cell:
			print(i)
		input('!!!!')
	return 1

def	get_synapses_pre_mtypes(filename):
	d_mtype_map, R_d_synapses, d_synapses, d_type_synapses = {}, {}, {}, {}

	with open(os.path.join(filename, 'mtype_map.tsv'), 'r') as fd:
		mtype_map = fd.read().split('\n')
		for line in mtype_map:
			line = line.split('\t')
			if len(line) > 1:
				d_mtype_map[int(line[0])] = line[1]

	'''
	d_mtype_map
	{0: 'L1_DAC', 1: 'L1_NGC-DA', 2: 'L1_NGC-SA', 3: 'L1_HAC', 4: 'L1_DLAC', 5: 'L1_SLAC', 6: 'L23_PC', 7: 'L23_MC', 8: 'L23_BTC', 9: 'L23_DBC', 10: 'L23_BP', 11: 'L23_NGC', 12:
	'''
	d_num_mtypes_synapses = {}
	with open(os.path.join(filename, 'synapses.tsv'), 'r') as fd:
		synapses = fd.read().split('\n')

		for line in synapses[1:]:

			line = line.split('\t')
			if len(line) > 1:
				synapse_id = int(line[0])
				pre_mtype = int(line[2])
				#If synapse_type < 100 the synapse is inhibitory, otherwise excitatory
				syn_pre_type = int(line[6])
				excitatory = 1
				if syn_pre_type < 100:
					excitatory = 0

				str_pre_mtype = d_mtype_map[pre_mtype]
				str_pre_mtype = str_pre_mtype.replace('-', '_')
				d_synapses[synapse_id] = str_pre_mtype
				
				R_d_synapses.setdefault(str_pre_mtype, [])
				R_d_synapses[str_pre_mtype].append(synapse_id)

				d_type_synapses[synapse_id] = excitatory

				d_num_mtypes_synapses.setdefault(str_pre_mtype, 0)
				d_num_mtypes_synapses[str_pre_mtype] += 1
	'''
	synapse_id = [0...100...N]

	##d_num_mtypes_synapses = {'L1_DAC':31, 'L4_SP':20}# mtype_pre:num | number of pre neuron mtypes
	## R_d_synapses = {'L1_DAC':0, 'L4_SP':1} mtype_pre:synapse_id
	##d_synapses = {0:'L1_DAC', 1:'L4_SP'} synapse_id:mtype_pre
	## d_type_synapses = {0:1, 1:0} synapse_id:excitatory | 0-inhibitory, 1-excitatory
	'''
	return d_num_mtypes_synapses, R_d_synapses, d_synapses, d_type_synapses


def get_lists():	
	str_id_n_id = {}
	gid_num_pre_syn_mtype, n_id_str_id, gid_syn_pre_cell_type, gid_R_syn_pre_cell_type = [], [], [], []
	
	ll = len(os.listdir('template'))
	for _asf, _template_paths in enumerate(os.listdir('template')):
		if _template_paths.endswith('.hoc'):
			post_nrn_str_id = _template_paths.split('.hoc')[0]
			print(ll-_asf, post_nrn_str_id)

			num_mtypes_synapses, R_mtype_synapses, mtype_synapses, ie_type_synapses = get_synapses_pre_mtypes(opj('synapses', post_nrn_str_id))

			## num_mtypes_synapses = {'L1_DAC':31, 'L4_SP':20}# mtype_pre:num | number of pre neuron mtypes
			## R_mtype_synapses = {'L1_DAC':0, 'L4_SP':1} mtype_pre:synapse_id
			## mtype_synapses = {0:'L1_DAC', 1:'L4_SP'} synapse_id:mtype_pre
			## ie_type_synapses = {0:1, 1:0} synapse_id:excitatory | 0-inhibitory, 1-excitatory

			str_id_n_id[post_nrn_str_id] = len(n_id_str_id) # {'L4_BP_bAC217_2': 0}
			n_id_str_id.append(post_nrn_str_id) # ['L4_BP_bAC217_2']
			gid_syn_pre_cell_type.append(mtype_synapses)
			gid_R_syn_pre_cell_type.append(R_mtype_synapses)
			gid_num_pre_syn_mtype.append(num_mtypes_synapses)

	return str_id_n_id, n_id_str_id, gid_syn_pre_cell_type, gid_R_syn_pre_cell_type, gid_num_pre_syn_mtype

def has_this_mtype(n_id_str_id, target_num, pre_syn_mtype, should_be_mtypes, pre_syn_cell_ids_available):
	tmp = []
	result = []
	for _target_num in range(len(n_id_str_id)): 
		if pre_syn_mtype in n_id_str_id[_target_num] and _target_num != target_num:
			tmp.append(_target_num)
	
	for i in range(should_be_mtypes):
		result.append(rnd.choice(tmp))
	return result

def has_this_mtype_z3(post_syn_cell_id, str_id_n_id, should_be_mtypes, pre_syn_cell_ids_available):

	if post_syn_cell_id in pre_syn_cell_ids_available:
		del pre_syn_cell_ids_available[post_syn_cell_id]
	result = []
	for i in range(should_be_mtypes):
		maxx = None
		pre_syn_cell_id = None
		for _pre_syn_cell_id, nums_left in pre_syn_cell_ids_available.items():
			if maxx is None or maxx < nums_left:
				maxx = nums_left
				pre_syn_cell_id = _pre_syn_cell_id
		if maxx is None:
			break
		del pre_syn_cell_ids_available[pre_syn_cell_id]
		result.append(str_id_n_id[pre_syn_cell_id])
	return result

def get_mtype(cell_id, inhib):
	for k, v in inhib.items():
		if cell_id.startswith(k):
			return k, v
	exit(f'no type for {cell_id}')

def duplicated_map_i_e(multiplier, folder_layers='AllLayers'):
	inhib, excit = get_inh_exi()
	print('------------in')
	print(inhib)
	print('------------ex')
	print(excit)
	print('------------')
	prior_etype = ["cNAC", "bNAC", "dNAC", "cAC", "bAC", "_"]
	cells = {'L1':[], 'L23':[],'L4':[],'L5':[],'L6':[]}
	for layer in os.listdir(folder_layers):
		print(layer)
		layer_cell = opj(folder_layers, layer)
		if os.path.isdir(layer_cell):
			for cell in os.listdir(layer_cell):
				if os.path.isdir(opj(layer_cell, cell)):
					cells[layer].append(cell)
					# print('	', cell)

	with open('/Users/nazargulya/Проекты/NEURO/NEURO_was/Connectome/div_resulting.json', 'r') as f:
		div = json.load(f)
	
	all_nrns = len(div)
	nums = {'L1':0, 'L23':0,'L4':0,'L5':0,'L6':0}
	inh_statistics = {'L1':0, 'L23':0,'L4':0,'L5':0,'L6':0}
	ex_statistics = {'L1':0, 'L23':0,'L4':0,'L5':0,'L6':0}

	for k1 in div.keys():
		k1 = k1.replace('-', '_')
		layer = k1.split('_')[0]
		this_inhib = -1
		for k2, v2 in inhib.items():
			if k1.find(k2) >= 0:
				this_inhib = v2
				break
		if this_inhib ==-1:
			exit(f'{k1} huy znaet')
		elif this_inhib == 0:
			ex_statistics[layer] += 1
		else:
			inh_statistics[layer] += 1
		nums[layer] += 1

	print('e', ex_statistics)
	print('i', inh_statistics)
	# weighted_ex = ex_statistics.copy()
	# for k, v in inh_statistics.items():
	# 	weighted_ex[k] = np.round(weighted_ex[k]/(v + weighted_ex[k]), 3)
	# 	print(k,weighted_ex[k])
	d_map = {}
	all_nums_ = 0
	for layer in inh_statistics.keys():
		d_map[layer] = {}
		print(layer)
		print(inh_statistics[layer], ex_statistics[layer])
		i_num = max([int(inh_statistics[layer] * multiplier), 1])
		e_num = max([int(ex_statistics[layer] * multiplier), 1])
		if layer == 'L1':
			e_num = 0
		print(i_num, e_num)
		
		i_tmp_etype = [0 for i in prior_etype]
		e_tmp_etype = [0 for i in prior_etype]
		prior_range = range(len(prior_etype))

		for available_cell in cells[layer]:
			this_inhib = -1
			available_cell_clean = available_cell.replace('-', '_')
			for c, i in inhib.items():
				if available_cell_clean.find(c)>=0:
					this_inhib = i
					break

			available_cell = available_cell.split('_')
			available_cell = '_'.join(available_cell[:-1])
			if this_inhib == 1:
				if i_tmp_etype[0] == 0:
					for i in prior_range:
						if i_tmp_etype[i] == 0 and available_cell.find(prior_etype[i]) >= 0:
							i_tmp_etype[i] = available_cell
							break
			elif this_inhib == 0 :
				if e_tmp_etype[0] == 0:
					for i in prior_range:
						if e_tmp_etype[i] == 0 and available_cell.find(prior_etype[i]) >= 0:
							e_tmp_etype[i] = available_cell
							break
			elif this_inhib == -1:
				exit(f'this_inhib = -1 {available_cell}')
			else:
				print('this_inhib', this_inhib)
				break
		for i in prior_range:
			if i_num > 0 and i_tmp_etype[i] != 0:
				all_nums_ += i_num
				d_map[layer][i_tmp_etype[i]] = max([i_num, 2])
				i_num = 0
			if e_num > 0 and e_tmp_etype[i] != 0:
				all_nums_ += e_num
				d_map[layer][e_tmp_etype[i]] = max([e_num, 2])
				e_num = 0
			

			print(prior_etype[i], '|', i_tmp_etype[i], e_tmp_etype[i])
		print('-----------------')
	print(d_map)
	print('all_nums_::', all_nums_)
	print('---------')
	# exit()
	return d_map

def get_more_genereal_i_e(inhib, excit, gid_num_pre_syn_mtype, gid_syn_pre_cell_type, gid_R_syn_pre_cell_type):
	'''
	gid_num_pre_syn_mtype:
	[0]{'L4_PC': 7, 'L4_SP': 2, 'L4_BTC': 3, 'L6_LBC': 2, 'L23_DBC': 3, 'L5_MC': 8, 'L5_NBC': 6}
	ie_gid_num_pre_syn_mtype:
	[0]{'L4_i': 317, 'L23_e': 231}

	gid_syn_pre_cell_type: as list
	[0]{0: 'L4_PC', 1: 'L4_PC', 2: 'L4_SP', 3: 'L4_SP', 4: 'L4_SP', 5: 'L4_SP', .....
	ie_gid_syn_pre_cell_type:
	[0]{0: 'L5_i', 1: 'L1_e', 2: 'L23_e', ....

	gid_R_syn_pre_cell_type:
	[gid]{'L4_PC': [0, 1, 6, 10, 11, 17, 25], 'L4_SP': [2, 3]}....
	ie_gid_R_syn_pre_cell_type:
	[0]{'L4_i': [0, 1, 6, 10, 11, 17, 25], 'L5_e': [2, 3]}
	'''
	ie_gid_num_pre_syn_mtype = []
	for d in gid_num_pre_syn_mtype:
		tmp_d = {}
		for pre_mtype, nums in d.items():
			name = '_e'
			if inhib[pre_mtype] == 1:
				name = '_i'
			name = pre_mtype.split('_')[0] + name
			tmp_d.setdefault(name, 0)
			tmp_d[name] += nums
		ie_gid_num_pre_syn_mtype.append(tmp_d)

	ie_gid_syn_pre_cell_type = []
	for d in gid_syn_pre_cell_type:
		tmp_d = {}
		for ix, pre_mtype in d.items():
			name = '_e'
			if inhib[pre_mtype] == 1:
				name = '_i'
			name = pre_mtype.split('_')[0] + name
			tmp_d[ix] = name
		ie_gid_syn_pre_cell_type.append(tmp_d)

	ie_gid_R_syn_pre_cell_type = []
	for d in gid_R_syn_pre_cell_type:
		tmp_d = {}
		for pre_mtype, list_ix_syn in d.items():
			name = '_e'
			if inhib[pre_mtype] == 1:
				name = '_i'
			name = pre_mtype.split('_')[0] + name
			tmp_d.setdefault(name, [])
			tmp_d[name].extend(list_ix_syn)
		ie_gid_R_syn_pre_cell_type.append(tmp_d)

	return ie_gid_num_pre_syn_mtype, ie_gid_syn_pre_cell_type, ie_gid_R_syn_pre_cell_type

def for_this_L_ie(pre_syn_L_ie, post_syn_num, n_id_str_id, inhib_gid):
	layer = pre_syn_L_ie.split('_')[0]
	i = 1 if pre_syn_L_ie.split('_')[-1] == 'i' else 0

	ans = []
	for ix, cell_id, in enumerate(n_id_str_id):
		if ix != post_syn_num:
			__layer = cell_id.split('_')[0]
			if __layer == layer and inhib_gid[ix] == i:
				ans.append(ix)
	return ans


def create_synapse_map_i_e(anatomy_path, folder=''):
	#difference num pre post
	inhib, excit = get_inh_exi()
	str_id_n_id, n_id_str_id, gid_syn_pre_cell_type, gid_R_syn_pre_cell_type, gid_num_pre_syn_mtype = get_lists()

	ie_gid_num_pre_syn_mtype, ie_gid_syn_pre_cell_type, ie_gid_R_syn_pre_cell_type = get_more_genereal_i_e(inhib, excit, \
		gid_num_pre_syn_mtype, gid_syn_pre_cell_type, gid_R_syn_pre_cell_type)
	
	inhib_gid = []
	for cell_id in n_id_str_id:
		ans = -1
		cell_id = cell_id.replace('-', '_')
		for cell, i in inhib.items():
			if cell_id.find(cell) >=0:
				ans = i
				break
		if ans == -1:
			exit(f'ans == -1, {cell_id}:cell_id')
		inhib_gid.append(ans)

	'''
	str_id_n_id:
	 {'L4_BP_bAC217_2': 0}
	n_id_str_id:
	 ['L4_BP_bAC217_2']
	

	gid_num_pre_syn_mtype:
	[0]{'L4_PC': 7, 'L4_SP': 2, 'L4_BTC': 3, 'L6_LBC': 2, 'L23_DBC': 3, 'L5_MC': 8, 'L5_NBC': 6}
	ie_gid_num_pre_syn_mtype:
	[0]{'L4_i': 317, 'L23_e': 231}

	gid_syn_pre_cell_type: as list
	[0]{0: 'L4_PC', 1: 'L4_PC', 2: 'L4_SP', 3: 'L4_SP', 4: 'L4_SP', 5: 'L4_SP', .....
	ie_gid_syn_pre_cell_type:
	[0]{0: 'L5_i', 1: 'L1_e', 2: 'L23_e', ....

	gid_R_syn_pre_cell_type:
	[gid]{'L4_PC': [0, 1, 6, 10, 11, 17, 25], 'L4_SP': [2, 3]}....
	ie_gid_R_syn_pre_cell_type:
	[0]{'L4_i': [0, 1, 6, 10, 11, 17, 25], 'L5_e': [2, 3]}
	'''
	synapses_map = {}
	num_cells = len(n_id_str_id)	

	for post_syn_num in range(num_cells):
		post_syn_cell_id = n_id_str_id[post_syn_num]
		print(num_cells-post_syn_num, post_syn_cell_id)

		pre_syn_mtype_num = ie_gid_num_pre_syn_mtype[post_syn_num].copy()#{'L4_i': 317, 'L23_e': 231}
		R_pre_syn_mtype = ie_gid_R_syn_pre_cell_type[post_syn_num].copy()#{'L4_i': [0, 1, 6, 10, 11, 17, 25], 'L5_e': [2, 3]}

		synapses_map[post_syn_cell_id] = gid_syn_pre_cell_type[post_syn_num].copy()#{0: 'L4_PC', 1: 'L4_PC', 2: 'L4_SP', 3: 'L4_SP', 4: 'L4_SP', 5: 'L4_SP', .....

		for pre_syn_L_ie, num_synapses in pre_syn_mtype_num.items():#{'L4_i': 317, 'L23_e': 231}
			
			available_gid = for_this_L_ie(pre_syn_L_ie, post_syn_num, n_id_str_id, inhib_gid)
			num_available = len(available_gid)
			
			print(pre_syn_L_ie, num_synapses, num_available)

			if num_available:
				drop = False
				if num_synapses <= num_available:
					drop = True

				for syn_id in R_pre_syn_mtype[pre_syn_L_ie]:
					cell_id = rnd.choice(available_gid)
					if drop:
						available_gid.remove(cell_id)
					synapses_map[post_syn_cell_id][syn_id] = n_id_str_id[cell_id]
			else:
				for syn_id in R_pre_syn_mtype[pre_syn_L_ie]:
					synapses_map[post_syn_cell_id][syn_id] = 'no-data'
					exit('no-data')
				# input('no-data')
	with open(f'result/synapses_map.json', 'w') as outfile:
		json.dump(synapses_map, outfile)

'''
#L1_SLAC_* no info in folder
	#L1_SAC_* no in anatomy, but have info /rename to SLAC_*

	L1_LAC_* no in anatomy /delete!

# L1_DLAC_cNAC no info in folder
# L1_DAC_cNAC have info /copy with new name DLAC!
'''


if __name__ == '__main__':
	should_mtypes = ['L1_DAC','L1_NGC_DA','L1_NGC_SA','L1_SLAC', 'L1_DLAC', 'L1_HAC','L23_PC','L23_MC','L23_BTC','L23_DBC','L23_BP','L23_NGC','L23_LBC','L23_NBC','L23_SBC','L23_ChC','L4_PC','L4_SP','L4_SS','L4_MC','L4_BTC','L4_DBC','L4_BP','L4_NGC','L4_LBC','L4_NBC','L4_SBC','L4_ChC','L5_TTPC1','L5_TTPC2','L5_UTPC','L5_STPC','L5_MC','L5_BTC','L5_DBC','L5_BP','L5_NGC','L5_LBC','L5_NBC','L5_SBC','L5_ChC','L6_TPC_L1','L6_TPC_L4','L6_UTPC','L6_IPC','L6_BPC','L6_MC','L6_BTC','L6_DBC','L6_BP','L6_NGC','L6_LBC','L6_NBC','L6_SBC','L6_ChC']

	opj = os.path.join
	multiplier = 0.1#426

	d_map = duplicated_map(multiplier=multiplier)
	load_duplicated_neurons_data(d_map, axon_use=False, circuit_folders='AllLayers')
	create_synapse_map_v2('BBPjson/pathways_anatomy_factsheets_simplified.json', 'BBPjson/')
	print("\n\nExecute this command:\nnrnivmodl mechanisms/")

	# if multiplier >= 0.01 or 1==1:
	# 	d_map = duplicated_map(multiplier=multiplier)
	# 	load_duplicated_neurons_data(d_map, axon_use=False, circuit_folders='AllLayers')
	# 	#### load_neurons_data(axon_use=False, circuit_folders=['AllLayers/L23'], name_cell='L23_SBC_dNAC222_3')
	# 	#### load_neurons_data(axon_use=False, circuit_folders=['MyLayers5'])	
	# 	create_synapse_map_v2('BBPjson/pathways_anatomy_factsheets_simplified.json', 'BBPjson/')
	# else:
	# 	d_map = duplicated_map_i_e(multiplier=multiplier)##all 90 12inh>>>>>>>>??????
	# 	load_duplicated_neurons_data(d_map, axon_use=False, circuit_folders='AllLayers')
	# 	create_synapse_map_i_e('BBPjson/pathways_anatomy_factsheets_simplified.json', 'BBPjson/')
