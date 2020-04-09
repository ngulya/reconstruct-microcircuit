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
				x.append(sec.x3d(i) + num_cell*1000)
				y.append(sec.y3d(i) + num_cell*1000)
				z.append(sec.z3d(i) + num_cell*1000)
				diam = sec.diam3d(i)
			# ax.plot3D(x, y, z, label=f'{sec}', linewidth=diam)
			ax.plot3D(x, y, z, color=color, linewidth=diam)
		for i in range(len(ids)):
			pos = additional[i]
			ax.scatter(pos[0], pos[1], pos[2], linewidth=4, label=ids[i])
	plt.legend()
	plt.show()

def get_synapses_pre_mtypes(filename):
	d_mtype_map, d_synapses = {}, {}

	with open(opj(filename, 'mtype_map.tsv'), 'r') as fd:
		mtype_map = fd.read().split('\n')
		for line in mtype_map:
			line = line.split('\t')
			if len(line) > 1:
				d_mtype_map[int(line[0])] = line[1]
	
	with open(opj(filename, 'synapses.tsv'), 'r') as fd:
		synapses = fd.read().split('\n')
		for line in synapses[1:]:
			line = line.split('\t')
			if len(line) > 1:
				synapse_id = int(line[0])
				pre_mtype = int(line[2])
				d_synapses[synapse_id] = d_mtype_map[pre_mtype]
	return d_synapses



def get_neurons_path(circuit_folder):
	l_neuron_folders = []
	for neuron_folders in os.listdir(circuit_folder):
		if os.path.isdir(opj(circuit_folder, neuron_folders)):
			l_neuron_folders.append(neuron_folders)
	return l_neuron_folders


def change_path_in_hoc(filename, text_to_search, text_to_replace):

	with fileinput.FileInput(filename, inplace=True) as file:
		for line in file:
			print(line.replace(text_to_search, text_to_replace), end='')
	return 1

def change_axons_settings(filename):
	first_replace = True
	with fileinput.FileInput(filename, inplace=True) as file:
		for line in file:
			if first_replace and line.find('replace_axon()') != -1:
				line = line.replace('replace_axon', '//replace_axon')
				first_replace = False
			print(line, end='')
	return 1

def load_neurons_data(file_loc, circuit_folder='MyLayers'):
	for folder in ['synapses', 'morphology', 'biophysics', 'template']:
		shutil.rmtree(folder, ignore_errors=True)
		os.makedirs(folder, exist_ok=True)
	shutil.rmtree('mechanisms', ignore_errors=True)
	
	have_mechanisms = False
	neuron_folders = get_neurons_path(circuit_folder)
	biophysic_names = {}
	for cell_id in neuron_folders:
		print(cell_id)
		folder = opj(circuit_folder, cell_id)
		with open(opj(folder, 'cellinfo.json'), 'r') as fd:
			cell_info = json.load(fd)
		biophysic_name = cell_info["cell name"].split('_')[0]

		biophysics_hoc_old = opj(folder, 'biophysics.hoc')
		biophysics_hoc_new = opj('biophysics', f'{cell_id}.hoc')
	
		if not biophysic_name in biophysic_names.keys():
			shutil.copy(biophysics_hoc_old, biophysics_hoc_new)
			biophysic_names[biophysic_name] = biophysics_hoc_new
		else:
			biophysics_hoc_new = biophysic_names[biophysic_name]


		morphology_old = opj(folder, 'morphology')
		morphology_new = opj('morphology', cell_id)
		shutil.copytree(morphology_old, morphology_new)

		morphology_hoc_old = opj(folder, 'morphology.hoc')
		morphology_hoc_new = opj(morphology_new, 'morphology.hoc')
		shutil.copy(morphology_hoc_old, morphology_hoc_new)
		change_path_in_hoc(opj(morphology_new, 'morphology.hoc'), 'morphology/', opj(file_loc, morphology_new+'/'))

		synapses_old = opj(folder, 'synapses')
		synapses_new = opj('synapses', cell_id)
		shutil.copytree(synapses_old, synapses_new)
		change_path_in_hoc(opj(synapses_new, 'synapses.hoc'), 'synapses/', opj(file_loc, synapses_new+'/'))
		change_path_in_hoc(opj(synapses_new, 'synapses.hoc'), 'synapse_list.append(synapse)', 'synapse.synapseID = synapse_id\nsynapse_list.append(synapse)')


		template_hoc_old = opj(folder, 'template.hoc')
		template_hoc_new = opj('template', f'{cell_id}.hoc')
		shutil.copy(template_hoc_old, template_hoc_new)
		change_path_in_hoc(template_hoc_new, 'morphology.hoc', opj(file_loc, morphology_hoc_new))
		change_path_in_hoc(template_hoc_new, 'biophysics.hoc', opj(file_loc, biophysics_hoc_new))
		change_path_in_hoc(template_hoc_new, 'synapses/synapses.hoc', opj(file_loc, synapses_new, 'synapses.hoc'))
		change_axons_settings(template_hoc_new)	

		if not have_mechanisms:
			shutil.copytree(opj(folder, 'mechanisms'), 'mechanisms')
			have_mechanisms = True
		else:
			diff_files = dircmp(opj(folder, 'mechanisms'), 'mechanisms').diff_files
			if len(diff_files) > 0:
				exit(f'Some files in {opj(folder, "mechanisms")} differences: {diff_files}')


	return 1

def get_lists():
	syn_pre_cell_type, template_cells, template_paths = [], [], []
	for _template_paths in os.listdir('template'):
		if _template_paths.endswith('.hoc'):
			cell_id = _template_paths.split('.hoc')[0]
			_template_paths = opj('template', _template_paths)
			with open(_template_paths, 'r') as fd:
				for line in fd.read().split('\n'):
					if line.find('begintemplate') != -1:
						template_cell = line.split('begintemplate')[1].strip()
						break

			syn_pre_cell_type.append(get_synapses_pre_mtypes(opj('synapses', cell_id)))
			template_paths.append(_template_paths)
			template_cells.append(template_cell)
	return template_cells, template_paths, syn_pre_cell_type

def load_neurons(template_cells, template_paths):
	neuron.h.load_file("stdrun.hoc")
	neuron.h.load_file("import3d.hoc")

	cells = []
	for num_cell in range(len(template_paths)):
		neuron.h.load_file(template_paths[num_cell])
		print(template_cells[num_cell])
		cells.append(getattr(neuron.h, template_cells[num_cell])(1))

	return cells

def run_vasya_run(tstop):
	print(f'\n\nRunning for {tstop} ms')
	tstart = time.time()
	leftsec = tstop
	num_windows = max([int(float(tstop)/50), 10])
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

def get_section_xyz_d(section, loc=0.5):
	loc = int(loc* section.n3d())
	return [section.x3d(loc), section.y3d(loc), section.z3d(loc), section.diam3d(loc)]


def connect_2_neurons(cells, syn_pre_cell_type):
	target = cells[0]#target
	source = cells[1]#source
	our_synapses = []
	for syn_num, syn_pre_type in syn_pre_cell_type[0].items():
		if syn_pre_type.startswith('L6'):
			our_synapses.append(syn_num)
	
	####### for sec in source.soma[0].wholetree():
	####### 	sec.nseg = 1
	####### for sec in target.soma[0].wholetree():
	####### 	sec.nseg = 1
	
	
	dist = []
	vol = []
	for i in range(0, 90, 10):
		dist.append(distance.euclidean(get_section_xyz_d(source.soma[0])[:3], get_section_xyz_d(source.axon[i])[:3]))
		vol.append(neuron.h.Vector().record(source.axon[i](0.5)._ref_v))
	
	iclamp = neuron.h.IClamp(source.soma[0](0.5))
	iclamp.delay = 100
	iclamp.dur = 300
	iclamp.amp = 0.9

	source_v = neuron.h.Vector().record(source.soma[0](0.5)._ref_v)
	source_axon_v = neuron.h.Vector().record(source.axon[10](0.5)._ref_v)

	# for i in source.axon:
	# 	print(i)
	# exit()
	target_v = neuron.h.Vector().record(target.soma[0](0.5)._ref_v)
	time_v = neuron.h.Vector().record(neuron.h._ref_t)
	# print(source.axon)
	# for i in source.axon:
	# 	print(source.axon[0](0.5))

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
	run_vasya_run(500.0)
	# plt.plot(time_v, source_v, label = 'source_v')
	for i in range(len(dist)):
		if i in [4,5]:
			plt.plot(time_v, vol[i], label = f'd:{dist[i]}')

	# plt.plot(time_v, source_axon_v, label = 'source_axon_v')
	# plt.plot(time_v, target_v, label = 'target_v')
	plt.legend()
	plt.show()

	return netcons


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
	file_loc = os.path.dirname(os.path.abspath(__file__))
	# load_neurons_data(file_loc)
	
	# os.system("nrnivmodl mechanisms")
	# exit()

	template_cells, template_paths, syn_pre_cell_type = get_lists()
	cells = load_neurons(template_cells, template_paths)
	# plotting(cells)

	print(template_cells)
	
	# # print(cells[0].soma[0](0.5).k_ion)
	# print('internal=', neuron.h.cai0_ca_ion, neuron.h.ki0_k_ion)
	# print('external=', neuron.h.cao0_ca_ion, neuron.h.ko0_k_ion)
	# exit()
	change_ion_conc(cells, syn_pre_cell_type)















