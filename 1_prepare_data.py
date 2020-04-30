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
rnd.seed(13)

def get_neurons_path(circuit_folder):
	l_neuron_folders = []
	for neuron_folders in os.listdir(circuit_folder):
		if os.path.isdir(opj(circuit_folder, neuron_folders)):
			l_neuron_folders.append(neuron_folders)
	return l_neuron_folders


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
		if nume > numm:
			e = rnd.choice(list(etype.keys()))
			if etype[e] > 1:
				etype[e] -= 1
		else:
			e = rnd.choice(list(etype.keys()))
			etype[e] += 1
		nume = sum(etype.values())
	return etype, mtype

def change_layer_cell_number(lst_layer_cell, multiplier, etype, mtype):
	
	#sum(etype.values()) != sum(mtype.values()) because of diffent microcircuits???
	multiplier_me = sum(mtype.values())/sum(etype.values())
	for _type in etype.keys():
		etype[_type] = int(multiplier_me * etype[_type])
	
	imtype, ietype = mtype.copy(), etype.copy()

	for m in mtype:
		mtype[m] = max([int(mtype[m]*multiplier), 1])
	for e in etype:
		etype[e] = max([int(etype[e]*multiplier), 1])

	cells = []

	for cell in lst_layer_cell:
		cell = cell.split('_')
		cell = '_'.join(cell[:-1])
		if not cell in cells:
			cells.append(cell)

	cell_mtype, cell_etype = [], []
	for metype in cells:
		for _mtype in mtype.keys():
			if _mtype in metype:
				_etype_with_shit = metype.replace(_mtype+'_', '')				
				for _etype in etype.keys():
					if _etype in _etype_with_shit:
						cell_etype.append(_etype)
						break
				cell_mtype.append(_mtype)
				break
	
	etype, mtype = get_equal_sum_mtype_etype(etype, mtype)
	MEtype, EMtype = init_weights(mtype, etype, cell_mtype, cell_etype)

	for m, e in MEtype.items():
		mnum = mtype[m]
		_s = sum(e.values())
		for _e, _num in e.items():
			v = max((int(_num/_s*mnum), 1))
			e[_e], EMtype[_e][m] = v, v
#		print(mtype[m],'	|',m, '	', e)#list(e.keys()), sum(e.values()))

	# for e, m in EMtype.items():
	# 	print(etype[e],'	|',e, '	', m)#ist(m.keys()), sum(m.values()))

	newmtype, newetype = get_new_statistics(mtype, etype, MEtype)
	print('---')
	print(imtype, 'imtype real')
	print(mtype, 'imtype real normalized')
	print(newmtype, 'imtype new')
	print('---')
	print(ietype, 'ietype real')
	print(etype, 'ietype real normalized')
	print(newetype, 'ietype new')
	return get_new_duplicated_list(MEtype, cells, cell_mtype, cell_etype)	

def get_cell(m, e, cells, cell_mtype, cell_etype):
	for i in range(len(cells)):
		if cell_mtype[i] == m and cell_etype[i] == e:
			return cells[i]

def get_new_duplicated_list(MEtype, cells, cell_mtype, cell_etype):
	resulting = {}
	for m, _e in MEtype.items():
		for e, nums in _e.items():
			cell_id = get_cell(m, e, cells, cell_mtype, cell_etype)

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


def duplicated_map(multiplier, json_folder='BBPjson', folder_layers='AllLayers'):
	'''
	multiplier = 0.1 means that only 10% from all 30k neurons will be modeleted
	'''
	cells = {'L1':[], 'L23':[],'L4':[],'L5':[],'L6':[]}
	for layer in os.listdir(folder_layers):
		layer_cell = opj(folder_layers, layer)
		if os.path.isdir(layer_cell):
			for cell in os.listdir(layer_cell):
				if os.path.isdir(opj(layer_cell, cell)):
					cells[layer].append(cell)
	
	with open(opj(json_folder, 'circuit_download.json'), 'r') as f:
		circuit_download = json.load(f)
	
	with open(opj(json_folder, 'layer_download.json'), 'r') as f:
		layer_download = json.load(f)

	num_cells, num_cells_in_crc = 0, 0
	# print(4873 +6081 +6453 +10080 + 144)
	d_map = {}
	for layer, lst_layer_cell in cells.items():
		# if layer != 'L23':
		# 	continue
		print(layer)
		_num_cells_in_crc = circuit_download["No. of neurons per layer"][layer]
		
		etype = layer_download[layer]["No. of neurons per electrical types"]
		mtype = layer_download[layer]["No. of neurons per morphological types"]
		print('_num_cells_in_crc:', _num_cells_in_crc, sum(mtype.values()), 'mtype')

		d_map[layer] = change_layer_cell_number(lst_layer_cell, multiplier, etype, mtype)
		print(d_map[layer])
		num_cells += sum(d_map[layer].values())
	print('will  be', num_cells)
	return d_map

def load_duplicated_neurons_data(d_map, axon_use=False, circuit_folders='AllLayers'):

	file_loc = os.path.dirname(os.path.abspath(__file__))
	
	for folder in ['synapses', 'morphology', 'biophysics', 'template', 'current_amps']:
		shutil.rmtree(folder, ignore_errors=True)
		os.makedirs(folder, exist_ok=True)
	shutil.rmtree('mechanisms', ignore_errors=True)
	
	have_mechanisms = False
	postfixs = ['_1','_2','_3','_4','_5']
	for layer in d_map.keys():
		for basic_cell_id, nums in d_map[layer].items():
			print(basic_cell_id, nums)
			for _num_cpy in range(nums):

				if _num_cpy < 5:
					pfx = postfixs[_num_cpy]
				else:
					pfx = rnd.choice(postfixs)
					
				folder = opj(circuit_folders, layer, basic_cell_id + pfx)
				with open(opj(folder, 'cellinfo.json'), 'r') as fd:
					cell_info = json.load(fd)

				cell_id = basic_cell_id.replace('-', '_') + pfx
				cell_id = f'{cell_id}_dplc_{_num_cpy}'
				# print('	', pfx, cell_id)
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
	return 1

def	get_synapses_pre_mtypes(filename):
	d_mtype_map, R_d_synapses, d_synapses, d_type_synapses = {}, {}, {}, {}

	with open(os.path.join(filename, 'mtype_map.tsv'), 'r') as fd:
		mtype_map = fd.read().split('\n')
		for line in mtype_map:
			line = line.split('\t')
			if len(line) > 1:
				d_mtype_map[int(line[0])] = line[1]
	
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
				
				str_mtype = d_mtype_map[pre_mtype]
				str_mtype = str_mtype.replace('-', '_')
				d_synapses[synapse_id] = str_mtype
				
				R_d_synapses.setdefault(str_mtype, [])
				R_d_synapses[str_mtype].append(synapse_id)

				d_type_synapses[synapse_id] = excitatory

				d_num_mtypes_synapses.setdefault(str_mtype, 0)
				d_num_mtypes_synapses[str_mtype] += 1
	'''
	synapse_id = [0...100...N]

	d_num_mtypes_synapses[m-type] = 2 	|nums pre_synaptic cell
	d_synapses[0..N] = 'L1_DAC'.. 		|m-type
	d_type_synapses[0..N] = 0/1  		|0-inhibitory, 1-excitatory
	'''
	return d_num_mtypes_synapses, R_d_synapses, d_synapses, d_type_synapses


def get_lists():	

	cid_gid = {}
	gid_num_pre_syn_mtype, gid_cid, gid_syn_pre_cell_type, gid_R_syn_pre_cell_type = [], [], [], []
	
	for _template_paths in os.listdir('template'):
		if _template_paths.endswith('.hoc'):
			cell_id = _template_paths.split('.hoc')[0]
			print(cell_id)

			num_mtypes_synapses, R_mtype_synapses, mtype_synapses, ie_type_synapses = get_synapses_pre_mtypes(opj('synapses', cell_id))
			'''
			L4_BP_bAC217_2
			cid_gid:
			 {'L4_BP_bAC217_2': 0}
			gid_cid:
			 ['L4_BP_bAC217_2']
			num_mtypes_synapses:
			 {'L4_PC': 7, 'L4_SP': 2, 'L4_BTC': 3, 'L6_LBC': 2, 'L23_DBC': 3, 'L5_MC': 8, 'L5_NBC': 6}
			mtype_synapses: as list
			 {0: 'L4_PC', 1: 'L4_PC', 2: 'L4_SP', 3: 'L4_SP', 4: 'L4_SP', 5: 'L4_SP', .....
			R_mtype_synapses:
			 {'L4_PC': [0, 1, 6, 10, 11, 17, 25], 'L4_SP': [2, 3]}....
			'''
			
			cid_gid[cell_id] = len(gid_cid)
			gid_cid.append(cell_id)
			gid_syn_pre_cell_type.append(mtype_synapses)
			gid_R_syn_pre_cell_type.append(R_mtype_synapses)
			gid_num_pre_syn_mtype.append(num_mtypes_synapses)

	return cid_gid, gid_cid, gid_syn_pre_cell_type, gid_R_syn_pre_cell_type, gid_num_pre_syn_mtype

def has_this_mtype(target_num, mtype, gid_cid):
	result = []
	for _target_num in range(len(gid_cid)): 
		if mtype in gid_cid[_target_num] and _target_num != target_num:
			result.append(_target_num)
			
	return result

def create_synapse_map():
	#difference num pre post
	cid_gid, gid_cid, gid_syn_pre_cell_type, gid_R_syn_pre_cell_type, gid_num_pre_syn_mtype = get_lists()

	'''
	cid_gid:
	 {'L4_BP_bAC217_2': 0}
	gid_cid:
	 ['L4_BP_bAC217_2']
	gid_num_pre_syn_mtype:
	[0]{'L4_PC': 7, 'L4_SP': 2, 'L4_BTC': 3, 'L6_LBC': 2, 'L23_DBC': 3, 'L5_MC': 8, 'L5_NBC': 6}
	gid_syn_pre_cell_type: as list
	[0]{0: 'L4_PC', 1: 'L4_PC', 2: 'L4_SP', 3: 'L4_SP', 4: 'L4_SP', 5: 'L4_SP', .....
	gid_R_syn_pre_cell_type:
	[gid]{'L4_PC': [0, 1, 6, 10, 11, 17, 25], 'L4_SP': [2, 3]}....
	'''
	synapses_map = {}
	for post_syn_num in range(len(gid_cid)):
		post_syn_cell_id = gid_cid[post_syn_num]
		pre_syn_mtype_num = gid_num_pre_syn_mtype[post_syn_num]#{'L4_PC': 7, 'L4_SP': 2, 'L4_BTC': 3, 'L6_LBC': 2, 'L23_DBC': 3, 'L5_MC': 8, 'L5_NBC': 6}
		
		pre_syn_mtype = gid_syn_pre_cell_type[post_syn_num]#{0: 'L4_PC', 1: 'L4_PC', 2: 'L4_SP', 3: 'L4_SP', 4: 'L4_SP', 5: 'L4_SP', ...
		R_pre_syn_mtype = gid_R_syn_pre_cell_type[post_syn_num]#{'L4_PC': [0, 1, 6, 10, 11, 17, 25], 'L4_SP': [2, 3]}....
		
		synapses_map[post_syn_cell_id] = pre_syn_mtype.copy()#{0: 'L4_PC', 1: 'L4_PC', 2: 'L4_SP', 3: 'L4_SP', 4: 'L4_SP', 5: 'L4_SP', .....

		for mtype, num_synapses in pre_syn_mtype_num.items():
			available_gid = has_this_mtype(post_syn_num, mtype, gid_cid)
			num_available = len(available_gid)
			
			if num_available:
				drop = False
				if num_synapses <= num_available:
					drop = True

				for syn_id in R_pre_syn_mtype[mtype]:
					cell_id = random.choice(available_gid)
					if drop:
						available_gid.remove(cell_id)
					synapses_map[post_syn_cell_id][syn_id] = gid_cid[cell_id]
			else:
				for syn_id in R_pre_syn_mtype[mtype]:
					synapses_map[post_syn_cell_id][syn_id] = 'no-data'

	with open(f'result/synapses_map.json', 'w') as outfile:
		json.dump(synapses_map, outfile)
	return synapses_map



if __name__ == '__main__':
	opj = os.path.join
	# d_map = duplicated_map(multiplier=0.005)
	# load_duplicated_neurons_data(d_map, axon_use=False, circuit_folders='AllLayers')
	
	# load_neurons_data(axon_use=False, circuit_folders=['AllLayers/L1'])
	# load_neurons_data(axon_use=False, circuit_folders=['MyLayers5'])
	# load_neurons_data(axon_use=False, circuit_folders=['MyLayers19'])
	# load_neurons_data(axon_use=False, circuit_folders=['AllLayers/L23'], name_cell='L23_SBC_dNAC222_3')
	
	# print("\n\nExecute this command:\nnrnivmodl mechanisms/")

	create_synapse_map()	

