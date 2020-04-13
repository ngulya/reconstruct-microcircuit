import os
import matplotlib
import numpy
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

def change_axons_settings(filename):
	first_replace = True
	with fileinput.FileInput(filename, inplace=True) as file:
		for line in file:
			if first_replace and line.find('replace_axon()') != -1:
				line = line.replace('replace_axon', '//replace_axon')
				first_replace = False
			print(line, end='')
	return 1

def load_neurons_data(axon_use=True, circuit_folders=['MyLayers']):
	d_nums = {'L1': 70, 'L23' :215, 'L4':230, 'L5':260, 'L6':260}
	d_nums_should = {'L1': 70, 'L23' :215, 'L4':230, 'L5':260, 'L6':260}

	file_loc = os.path.dirname(os.path.abspath(__file__))
	
	for folder in ['synapses', 'morphology', 'biophysics', 'template']:
		shutil.rmtree(folder, ignore_errors=True)
		os.makedirs(folder, exist_ok=True)
	shutil.rmtree('mechanisms', ignore_errors=True)
	
	have_mechanisms = False
	neuron_parent_folders = []
	for folder in circuit_folders:
		neuron_parent_folders.append(get_neurons_path(folder))

	for _num_fold, neuron_folders in enumerate(neuron_parent_folders):
		for cell_id in neuron_folders:
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
				change_axons_settings(template_hoc_new)

			if not have_mechanisms:
				shutil.copytree(opj(folder, 'mechanisms'), 'mechanisms')
				have_mechanisms = True
			else:
				new_files = dircmp(opj(folder, 'mechanisms'), 'mechanisms').left_only
				for file in new_files:
					print('\n\n\n\n\n\n\nnrnivmodl mechanisms')
					shutil.copy(opj(folder, 'mechanisms', file), opj('mechanisms', file))

	return 1

def reduce_layer_cell_number(lst_layer_cell, will_be_num, etype, mtype):

	multiplier_me = sum(mtype.values())/sum(etype.values())

	for _type in etype.keys():
		asad = etype[_type]
		etype[_type] = int(multiplier_me * etype[_type])#sum(etype.values()) != sum(mtype.values()) because of diffent microcircuits?

	cells = []
	for cell in lst_layer_cell:
		cell = cell.split('_')
		cell = '_'.join(cell[:-1])
		if not cell in cells:
			cells.append(cell)
	
	cell_mtype = []
	cell_etype = []
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

	etype_our, mtype_our = etype.copy(), mtype.copy()
	for k in etype_our.keys():
		etype_our[k] = 0
	for k in mtype_our.keys():
		mtype_our[k] = 0
	
	metype = {}
	for i in range(len(cells)):
		_etype = cell_etype[i]
		_mtype = cell_mtype[i]

		etype_our[_etype] += 1
		mtype_our[_mtype] += 1
		metype.setdefault(_mtype, [])
		metype[_mtype].append(_etype)

		# print(f'{cells[i]}	| {cell_mtype[i]} |	{cell_etype[i]}')
	# print(len(cells), len(numpy.unique(cell_mtype)), len(numpy.unique(cell_etype)))
	
	for k, v in metype.items():
		print(k, v)

	res_metype = {}
	for i in range(len(cells)):
		_mtype = cell_mtype[i]
		_etype = cell_etype[i]
		multiplier_m = mtype[_mtype]/len(metype[_mtype])
		print(f'	{_mtype}:should be {mtype[_mtype]}  but we have  {len(metype[_mtype])} | {multiplier_m}')
		
		res_metype.setdefault(_mtype, [])
		for _i in range(int(multiplier_m)):
			res_metype[_mtype].append(_etype)


	new_mtype_our, new_etype_our = mtype_our.copy(), etype_our.copy()
	for k in new_etype_our.keys():
		new_etype_our[k] = 0
	for k in new_mtype_our.keys():
		new_mtype_our[k] = 0
	
	for _mtype, _etype in res_metype.items():
		new_mtype_our[_mtype] = len(_etype)
		for _e in _etype:
			new_etype_our[_e] += 1

	print(new_mtype_our, 'New')
	print(mtype_our, 'Our')
	print(mtype)
	print()
	print(new_etype_our, 'New')
	print(etype_our, 'Our')
	print(etype)

	exit()


def duplicated_map(multiplier=0.1, json_folder='BBPjson', folder_layers='AllLayers'):
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
	for layer, lst_layer_cell in cells.items():
		_num_cells = len(lst_layer_cell)
		_num_cells_in_crc = circuit_download["No. of neurons per layer"][layer]
		# print(layer)
		will_be_num = _num_cells_in_crc*0.1
		etype = layer_download[layer]["No. of neurons per electrical types"]
		mtype = layer_download[layer]["No. of neurons per morphological types"]
		reduce_layer_cell_number(lst_layer_cell, will_be_num, etype, mtype)
		# if will_be_num <= _num_cells:
		# 	etype = layer_download[layer]["No. of neurons per electrical types"]
		# 	mtype = layer_download[layer]["No. of neurons per morphological types"]
		# 	reduce_layer_cell_number(lst_layer_cell, will_be_num, etype, mtype)
		# 	exit()
		num_cells += _num_cells
		num_cells_in_crc += _num_cells_in_crc
		print()
	print(num_cells, num_cells_in_crc)

	exit()

if __name__ == '__main__':
	opj = os.path.join
	duplicated_map()
	exit()
	load_neurons_data(axon_use=False, circuit_folders=['AllLayers/L1', 'AllLayers/L23', 'AllLayers/L4', 'AllLayers/L5', 'AllLayers/L6'])
	# load_neurons_data(axon_use=False, circuit_folders=['AllLayers/L1'])	
	# load_neurons_data(axon_use=False, circuit_folders=['MyLayers19'])
	print("\n\nExecute this command:\nnrnivmodl mechanisms/")

	

