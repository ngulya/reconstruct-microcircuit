import os
import neuron
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

# from neuron.units import ms, mV
# import neurom
# import neurom.viewer
# neurom.viewer.draw(neurom.load_neuron('morphology/sm090918b1-3_idA_-_Scale_x1.000_y1.025_z1.000.asc'), mode='3d')
# plt.show()
#############clear; nrniv -pyexe pyexe run.py 

def plotting(cells, additional=[], ids=[]):
	###change cooord
	# for num, c in enumerate(cells):
	# 	for sec in c.soma[0].wholetree():
	# 		for i in range(sec.n3d()):
	# 			sec.pt3dchange(i,num*1000 + sec.x3d(i),num*1000 + sec.y3d(i),num*1000 + sec.z3d(i),sec.diam3d(i))
	
	fig = plt.figure()
	ax = plt.axes(projection='3d')
	for num_cell, cell in enumerate(cells):
		all_tree = cell.soma[0].wholetree()
		for sec in all_tree:
			print('\n', sec)
			if 'dend' in str(sec):
				color='red'
			elif 'axon' in str(sec):
				color = 'blue'
			else:
				color='black'
			
			x,y,z=[],[],[]
			diam = 0
			for i in range(sec.n3d()):
				print(sec.x3d(i), sec.y3d(i), sec.z3d(i), sec.diam3d(i))
				x.append(sec.x3d(i) + num_cell*1)
				y.append(sec.y3d(i) + num_cell*1)
				z.append(sec.z3d(i) + num_cell*1)
				diam = sec.diam3d(i)
			# ax.plot3D(x, y, z, label=f'{sec}', linewidth=diam)
			ax.plot3D(x, y, z, color=color, linewidth=diam)
		for i in range(len(ids)):
			pos = additional[i]
			ax.scatter(pos[0], pos[1], pos[2], linewidth=4, label=ids[i])
	plt.legend()
	plt.show()

def get_synapses_pre_mtypes(filename):
	d_mtype_map, d_synapses, d_type_synapses = {}, {}, {}

	with open(opj(filename, 'mtype_map.tsv'), 'r') as fd:
		mtype_map = fd.read().split('\n')
		for line in mtype_map:
			line = line.split('\t')
			if len(line) > 1:
				d_mtype_map[int(line[0])] = line[1]
	
	d_num_mtypes_synapses = {}
	with open(opj(filename, 'synapses.tsv'), 'r') as fd:
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
				d_synapses[synapse_id] = str_mtype
				d_type_synapses[synapse_id] = excitatory
				d_num_mtypes_synapses.setdefault(str_mtype, 0)
				d_num_mtypes_synapses[str_mtype] += 1

	return d_num_mtypes_synapses, d_synapses, d_type_synapses



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

def change_axons_settings(filename):
	first_replace = True
	with fileinput.FileInput(filename, inplace=True) as file:
		for line in file:
			if first_replace and line.find('replace_axon()') != -1:
				line = line.replace('replace_axon', '//replace_axon')
				first_replace = False
			print(line, end='')
	return 1

def load_neurons_data(circuit_folders=['MyLayers']):
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
			change_axons_settings(template_hoc_new)

			if not have_mechanisms:
				shutil.copytree(opj(folder, 'mechanisms'), 'mechanisms')
				have_mechanisms = True
			else:
				new_files = dircmp(opj(folder, 'mechanisms'), 'mechanisms').left_only
				for file in new_files:
					print('\n\n\n\n\n\n\nnrnivmodl mechanisms')
					shutil.copy(opj(folder, 'mechanisms', file), opj('mechanisms', file))

				# if len(diff_files) > 0:
				# 	exit(f'Some files in {opj(folder, "mechanisms")} differences: {diff_files}')


	return 1

def get_lists():
	dict_with_info = {}
	syn_pre_cell_type, template_cell_ids = [], []
	for _template_paths in os.listdir('template'):
		if _template_paths.endswith('.hoc'):
			cell_id = _template_paths.split('.hoc')[0]
			dict_with_info[cell_id] = {'num_synapses':0, 'mtype_synapses':{}, 'ie_type_synapses':0}
			
			num_mtypes_synapses, mtype_synapses, ie_type_synapses = get_synapses_pre_mtypes(opj('synapses', cell_id))
			dict_with_info[cell_id]['num_synapses'] = len(mtype_synapses)
			dict_with_info[cell_id]['num_mtypes_synapses'] = num_mtypes_synapses
			dict_with_info[cell_id]['mtype_synapses'] = mtype_synapses
			dict_with_info[cell_id]['ie_type_synapses'] = ie_type_synapses
			syn_pre_cell_type.append(mtype_synapses)
			template_cell_ids.append(cell_id)

	return template_cell_ids, syn_pre_cell_type, dict_with_info

def load_neurons(template_cell_ids):
	neuron.h.load_file("stdrun.hoc")
	neuron.h.load_file("import3d.hoc")

	cells = []
	alls = len(template_cell_ids)
	for cell_id in template_cell_ids:
	
		neuron.h.load_file(opj('template', f'{cell_id}.hoc'))
		print(cell_id, 'left:', alls)
		cells.append(getattr(neuron.h, cell_id)(1))
		alls -= 1
	return cells

def run_vasya_run(tstop):
	tstop = 100
	print(f'\n\nRunning for {tstop} ms')
	tstart = time.time()
	leftsec = tstop
	num_windows = max([int(float(tstop)/50), 100])
	windows = numpy.linspace(0,tstop,num_windows)

	for iround, neuron.h.tstop in enumerate(windows[1:]):
		print(f'\rLeft {int(tstop-neuron.h.tstop)} ms in simulation || {leftsec} sec', end='    ')
		if iround == 0:
			neuron.h.run()
		else:
			neuron.h.continuerun(neuron.h.tstop)
			tspent = time.time() - tstart
			leftsec = numpy.round(tspent*(tstop-neuron.h.tstop)/neuron.h.tstop, 1)

	print(f'\rElapsed time {numpy.round(time.time() - tstart, 2)} sec .......')
	# exit()
def get_section_xyz_d(section, loc=0.5):
	loc = int(loc* section.n3d())
	return [section.x3d(loc), section.y3d(loc), section.z3d(loc), section.diam3d(loc)]


def connect_2_neurons(cells, syn_pre_cell_type):
	target = cells[0]#target
	source = cells[1]#source
	our_synapses = []
	for syn_num, syn_pre_type in syn_pre_cell_type[0].items():
		if syn_pre_type.startswith('L1'):
			our_synapses.append(syn_num)
	
	####### for sec in source.soma[0].wholetree():
	####### 	sec.nseg = 1
	####### for sec in target.soma[0].wholetree():
	####### 	sec.nseg = 1
	
	
	dist = []
	vol = []
	for i in range(0, 50, 10):
		dist.append(distance.euclidean(get_section_xyz_d(source.soma[0])[:3], get_section_xyz_d(source.axon[i])[:3]))
		vol.append(neuron.h.Vector().record(source.axon[i](0.5)._ref_v))
	
	iclamp = neuron.h.IClamp(source.soma[0](0.5))
	iclamp.delay = 100
	iclamp.dur = 300
	iclamp.amp = 0.1

	source_v = neuron.h.Vector().record(source.soma[0](0.5)._ref_v)
	source_axon_v = neuron.h.Vector().record(source.axon[10](0.5)._ref_v)

	target_v = neuron.h.Vector().record(target.soma[0](0.5)._ref_v)
	time_v = neuron.h.Vector().record(neuron.h._ref_t)

	# netcons = []
	# for syn in target.synapses.synapse_list:
	# 	if syn.synapseID in our_synapses:
	# 		if not 'ProbAMPANMDA_EMS' in str(syn):
	# 			continue
	# 		print(syn)
	# 		nc = neuron.h.NetCon(source.axon[10](0.5)._ref_v, syn, sec=source.axon[10])
	# 		nc.weight[0] = 0.05
	# 		nc.delay = 1
	# 		netcons.append(nc)
	run_vasya_run(100.0)
	# plt.plot(time_v, source_v, label = 'source_v')
	for i in range(len(dist)):
		if i in [4,5]:
			plt.plot(time_v, vol[i], label = f'd:{dist[i]}')

	# plt.plot(time_v, source_axon_v, label = 'source_axon_v')
	# plt.plot(time_v, target_v, label = 'target_v')
	plt.legend()
	plt.show()

	return 1


def change_ion_conc(cells, syn_pre_cell_type):
	source = cells[0]

	print('[ca] inter:', neuron.h.cai0_ca_ion, '	', neuron.h.cao0_ca_ion, 'ext')
	print('soma [ca] inter:', source.soma[0](0.5).cai, '	', source.soma[0](0.5).cao, 'ext')
	print(source.soma[0](0.5).eca, 'eca')
	neuron.h.cai0_ca_ion = 5e-05
	neuron.h.cao0_ca_ion = 2
	neuron.h.stdinit()

	print('[ca] inter:', neuron.h.cai0_ca_ion, '	', neuron.h.cao0_ca_ion, 'ext')
	print('soma [ca] inter:', source.soma[0](0.5).cai, '	', source.soma[0](0.5).cao, 'ext')
	print(source.soma[0](0.5).eca, 'eca')

	iclamp = neuron.h.IClamp(source.dend[1](0.5))
	iclamp.delay = 100
	iclamp.dur = 500
	iclamp.amp = 100e-3


	source_v = neuron.h.Vector().record(source.soma[0](0.5)._ref_v)
	source_dend_v = neuron.h.Vector().record(source.dend[1](0.5)._ref_v)
	source_axon_v = neuron.h.Vector().record(source.axon[0](0.5)._ref_v)
	time_v = neuron.h.Vector().record(neuron.h._ref_t)

	run_vasya_run(700.0)

	# plt.plot(time_v, source_axon_v, label = 'source_axon_v')
	plt.plot(time_v, source_v, label = 'source_v')
	plt.plot(time_v, source_dend_v, label = 'source_dend_v')
	
	
	plt.legend()
	plt.show()

	return 1

if __name__ == '__main__':
	opj = os.path.join
	
	# load_neurons_data(circuit_folders=['AllLayers/L1', 'AllLayers/L23', 'AllLayers/L4', 'AllLayers/L5', 'AllLayers/L6'])
	# load_neurons_data(circuit_folders=['AllLayers/L1'])
	load_neurons_data()
	# # # os.system("nrnivmodl mechanisms")
	# exit('=')

	template_cell_ids, syn_pre_cell_type, dict_with_info = get_lists()
	print(template_cell_ids)

	# with open('synapse_info.json', 'w') as outfile:
	# 	json.dump(dict_with_info, outfile)
	# exit()
	cells = load_neurons(template_cell_ids)
	plotting(cells)
	
	# for i in range(len(template_paths)):
	# 	print(template_paths[i])
	# 	print('	', syn_pre_cell_type[i])
	# 	print()
	# exit()
	
	# change_ion_conc(cells, syn_pre_cell_type)
	# connect_2_neurons(cells, syn_pre_cell_type)















