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
import random as rnd


rnd.seed(13)

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

def get_lists(save_dict):
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

	if save_dict:
		with open('synapse_info.json', 'w') as outfile:
			json.dump(dict_with_info, outfile)
	return template_cell_ids, syn_pre_cell_type, dict_with_info

def load_cells(template_cell_ids):
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

def run_run_run(tstop):
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
	# exit()
def get_section_xyz_d(section, loc=0.5):
	loc = int(loc* section.n3d())
	return [section.x3d(loc), section.y3d(loc), section.z3d(loc), section.diam3d(loc)]

def connect_neurons(cells, template_cell_ids, syn_pre_cell_type, axon_use=True):
	#[L4_BP_bAC217_2[0], L5_BP_bAC217_3[0], L1_DAC_cNAC187_3[0], L23_BP_bNAC219_3[0], L6_BP_bAC217_3[0]]
	netcons = []
	def get_short_inf(inf, s_list):
		d = {'L1':0, 'L23':0,'L4':0,'L5':0,'L6':0}
		for syn in s_list:
			l = inf[int(syn.synapseID)].split('_')[0]
			d[l] += 1
		return d

	for num_t, target in enumerate(cells):
		_inf = get_short_inf(syn_pre_cell_type[num_t], target.synapses.synapse_list)
		print('postsyn:', target, '	num pre syn:', len(target.synapses.synapse_list), _inf)
		for num_s, source in enumerate(cells):
			if target == source:
				continue
			print('	presyn:', source, end = '->	')
			max_axon = len(source.axon) - 1
			layer_source = template_cell_ids[num_s].split('_')[0]
			info_for_target = syn_pre_cell_type[num_t]
			_cnt = 0
			
			for syn in target.synapses.synapse_list:
				could_be_pre_type = info_for_target[int(syn.synapseID)]
				pre_layer = could_be_pre_type.split('_')[0]
				if pre_layer == layer_source:
					print(could_be_pre_type, end=' ')
					if _cnt > 15:
						_cnt = 0
						print('\n			', end='')

					source_location = source.soma[0]
					if axon_use:
						axon_pos = rnd.randint(0, max_axon)
						source_location = source.axon[axon_pos]

					nc = neuron.h.NetCon(source_location(0.5)._ref_v, syn, sec=source_location)
					nc.weight[0] = 0.05
					nc.delay = 1
					netcons.append(nc)
					_cnt += 1
			print()
		print('all presyn connected to postsyn:', target, end='\n\n')
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

	run_run_run(700.0)

	# plt.plot(time_v, source_axon_v, label = 'source_axon_v')
	plt.plot(time_v, source_v, label = 'source_v')
	plt.plot(time_v, source_dend_v, label = 'source_dend_v')
	
	
	plt.legend()
	plt.show()

	return 1

def check_axon_voltage(cells):
	cell = cells[0]
	iclamp = neuron.h.IClamp(cell.soma[0](0.5))
	iclamp.delay = 100
	iclamp.dur = 500
	iclamp.amp = 100e-3
	source_v = neuron.h.Vector().record(cell.soma[0](0.5)._ref_v)
	lv = []
	step = 80
	for i in range(10, len(cell.axon), step):
		lv.append(neuron.h.Vector().record(cell.axon[i](0.5)._ref_v))
	time_v = neuron.h.Vector().record(neuron.h._ref_t)
	run_run_run(700.0)

	for i, __ in enumerate(range(10, len(cell.axon), step)):
		plt.plot(time_v, lv[i], label = f'axon[{i}]')
	plt.plot(time_v, source_v, label = 'source_v')
	plt.legend()
	plt.show()

def set_stimuls_and_go(config, cells, template_cell_ids, netcons):
	# iclamp = neuron.h.IClamp(cells[3].soma[0](0.5))
	# iclamp.delay = 100
	# iclamp.dur = 300
	# iclamp.amp = 0.1

	iclamps = []
	for cell in cells:
		iclamp = neuron.h.IClamp(cell.soma[0](0.5))
		iclamp.delay = config['iclamp_delay']
		iclamp.dur = config['iclamp_dur']
		iclamp.amp = config['iclamp_amp']
		iclamps.append(iclamp)

	list_v = []
	for cell in cells:
		list_v.append(neuron.h.Vector().record(cell.soma[0](0.5)._ref_v))

	time_v = neuron.h.Vector().record(neuron.h._ref_t)
	run_run_run(config['tstop'])

	d_result = {'t':list(time_v), 'template_cell_ids':template_cell_ids, 'record':{}}

	for num_cell in range(len(list_v)):
		d_result['record'][template_cell_ids[num_cell]] = list(list_v[num_cell])
	
	with open(f'result/single.json', 'w') as outfile:
		json.dump(d_result, outfile)
	print(f'result/single.json')
	
	# for i in range(len(cells)):
	# 	label = str(cells[i])
	# 	if i == 3:
	# 		label = 'Main: ' + label
	# 	plt.plot(time_v, list_v[i], label = label)
	# plt.legend()
	# plt.show()
	

if __name__ == '__main__':

	opj = os.path.join
	config = {'tstop':600, 'iclamp_delay':50, 'iclamp_dur':100, 'iclamp_amp':0.1}
	
	template_cell_ids, syn_pre_cell_type, dict_with_info = get_lists(save_dict=False)
	print(template_cell_ids)
	cells = load_cells(template_cell_ids)
	
	# plotting(cells)
	# check_axon_voltage(cells)
	# change_ion_conc(cells, syn_pre_cell_type)
	netcons = connect_neurons(cells, template_cell_ids, syn_pre_cell_type, False)
	set_stimuls_and_go(config, cells, template_cell_ids, netcons)
