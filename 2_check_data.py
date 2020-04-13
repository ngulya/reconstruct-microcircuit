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
import random


random.seed(13)

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
	d_mtype_map, R_d_synapses, d_synapses, d_type_synapses = {}, {}, {}, {}

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
				
				R_d_synapses.setdefault(str_mtype, [])
				R_d_synapses[str_mtype].append(synapse_id)

	'''
	synapse_id = [0...100...N]
	d_num_mtypes_synapses[m-type] = 2 	|nums pre_synaptic cell
	R_d_synapses['L1_DAC'] = [0..N] 	|m-type
	d_synapses[0..N] = 'L1_DAC'.. 		|m-type
	d_type_synapses[0..N] = 0/1  		|0-inhibitory, 1-excitatory
	'''
	return d_num_mtypes_synapses, R_d_synapses, d_synapses, d_type_synapses


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

def load_cells(template_cell_ids):
	neuron.h.load_file("stdrun.hoc")
	neuron.h.load_file("import3d.hoc")

	cells = []
	cells_map = {}
	alls = len(template_cell_ids)
	for cell_id in template_cell_ids:
	
		neuron.h.load_file(opj('template', f'{cell_id}.hoc'))
		print(cell_id, 'left:', alls)
		_c = getattr(neuron.h, cell_id)(1)
		cells.append(_c)
		cells_map[cell_id] = _c
		alls -= 1
	return cells, cells_map

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


def get_lists(save_dict):
	dict_with_info = {}
	syn_pre_cell_type, template_cell_ids, R_syn_pre_cell_type, num_pre_syn_type = [], [], [], []
	for _template_paths in os.listdir('template'):
		if _template_paths.endswith('.hoc'):
			cell_id = _template_paths.split('.hoc')[0]
			dict_with_info[cell_id] = {'num_synapses':0, 'mtype_synapses':{}, 'ie_type_synapses':0}
			
			num_mtypes_synapses, R_mtype_synapses, mtype_synapses, ie_type_synapses = get_synapses_pre_mtypes(opj('synapses', cell_id))
			'''
			synapse_id = [0...100...N]

			num_mtypes_synapses[m-type] = 2 				|(dict) nums pre_synaptic cell
			R_mtype_synapses['L1_DAC'] = [0,1,10]			|(dict) m-type
			mtype_synapses[0...100...N] = 'L1_DAC'.. 		|(list) m-type
			ie_type_synapses[0...100...N] = 0/1  			|(list) 0-inhibitory, 1-excitatory
			'''
			dict_with_info[cell_id]['num_synapses'] = len(mtype_synapses)
			dict_with_info[cell_id]['num_mtypes_synapses'] = num_mtypes_synapses
			dict_with_info[cell_id]['mtype_synapses'] = mtype_synapses
			dict_with_info[cell_id]['ie_type_synapses'] = ie_type_synapses
			
			R_syn_pre_cell_type.append(R_mtype_synapses)
			syn_pre_cell_type.append(mtype_synapses)
			num_pre_syn_type.append(num_mtypes_synapses)
			template_cell_ids.append(cell_id)

	if save_dict:
		with open('synapse_info.json', 'w') as outfile:
			json.dump(dict_with_info, outfile)

	return template_cell_ids, \
			R_syn_pre_cell_type, \
			syn_pre_cell_type, \
			num_pre_syn_type, \
			dict_with_info

def _has_this_mtype(gid_cid, target_num, mtype):
	result = []
	for _target_num in range(len(gid_cid)): 
		if mtype in gid_cid[_target_num] and _target_num != target_num:
			result.append(_target_num)			
	return result

def _create_synapses(gid_cid, gid_num_pre_syn_mtype, gid_syn_pre_cell_type, gid_R_syn_pre_cell_type):
	# gid_cid[0]= 'L1_HAC_bNAC219_3'
	# gid_num_pre_syn_mtype[0]= {'L1_HAC': 19, 'L23_NBC': 32, 'L23_PC': 21, 'L23_MC': 69,'L6_NBC_TP': 3}
	# gid_syn_pre_cell_type[0]{0:'L1_DAC'}.. 		|m-type
	# gid_R_syn_pre_cell_type[0]{'L1_DAC': [0,1,3,4]}
	
	synapses_map = {}
	for target_num in range(len(gid_cid)):
		target_cell_id = gid_cid[target_num]
		pre_syn_mtype = gid_num_pre_syn_mtype[target_num]
		mtype_synapses = gid_syn_pre_cell_type[target_num]
		R_mtype_synapses = gid_R_syn_pre_cell_type[target_num]
		
		synapses_map[target_cell_id] = mtype_synapses.copy()#{0:'L1_DAE'}

		for mtype, num_synapses in pre_syn_mtype.copy().items():
			available = _has_this_mtype(gid_cid, target_num, mtype)
			num_available = len(available)
			
			if num_available:
				drop = False
				if num_synapses <= num_available:
					drop = True
				
				for syn_id in R_mtype_synapses[mtype]:
					cell_id = random.choice(available)
					if drop:
						available.remove(cell_id)
					synapses_map[target_cell_id][syn_id] = gid_cid[cell_id]
			else:
				for syn_id in R_mtype_synapses[mtype]:
					synapses_map[target_cell_id][syn_id] = 'no-data'

		# new_R_mtype_synapses = R_mtype_synapses.copy()
		# for pre_syn_mtype in R_mtype_synapses:
		# 	res = []
		# 	for syn_id in R_mtype_synapses[pre_syn_mtype]:
		# 		res.append(synapses_map[target_cell_id][syn_id])
		# 	l = list(numpy.unique(res))
		# 	new_R_mtype_synapses[pre_syn_mtype] = [l, len(res)]
		# asgid_cid = gid_cid.copy()
		# asgid_cid.sort()
		# l = 'L1'
		# for i in asgid_cid:
		# 	ll = i.split('_')[0]
		# 	if ll == l:
		# 		print('  ', i, end=' ')
		# 	else:
		# 		l = ll
		# 		print('\n  ', i, end=' ')
		# print()
		# # print(asgid_cid)
		# print('target_cell_id=', target_cell_id)
		# print('--')
		# for k, v in new_R_mtype_synapses.items():
		# 	print('	', k, v)		
		# input('------')

	return synapses_map


def get_short_inf(inf, s_list):
	d = {'L1':0, 'L23':0,'L4':0,'L5':0,'L6':0}
	for syn in s_list:
		l = inf[int(syn.synapseID)].split('_')[0]
		d[l] += 1
	return d

def connect_neurons(synapses_map, cells_map, template_cell_ids, syn_pre_cell_type, axon_use=True):
	print(synapses_map)
	netcons = []

	for target_cell_id, connect_map in synapses_map.items():
		post_syn = cells_map[target_cell_id]
		
		index_target = None
		for _i in range(len(template_cell_ids)):
			if template_cell_ids[_i] == target_cell_id:
				index_target = _i
				break

		_inf = get_short_inf(syn_pre_cell_type[index_target], post_syn.synapses.synapse_list)
		print('post_syn:', post_syn, '	num pre_syn:', len(post_syn.synapses.synapse_list))
		con, nocon = {'L1':0, 'L23':0,'L4':0,'L5':0,'L6':0}, {'L1':0, 'L23':0,'L4':0,'L5':0,'L6':0}
		for syn in post_syn.synapses.synapse_list:
			pre_syn_cell_mtype = connect_map[int(syn.synapseID)]

			pre_mtype_layer = syn_pre_cell_type[index_target][int(syn.synapseID)].split('_')[0]
			if pre_syn_cell_mtype == 'no-data':
				nocon.setdefault(pre_mtype_layer, 0)
				nocon[pre_mtype_layer] += 1
				continue

			con.setdefault(pre_mtype_layer, 0)
			con[pre_mtype_layer] += 1
			
			pre_syn = cells_map[pre_syn_cell_mtype]
			
			nc = neuron.h.NetCon(pre_syn.soma[0](0.5)._ref_v, syn, sec=pre_syn.soma[0])
			nc.weight[0] = 0.05
			nc.delay = 1
			netcons.append(nc)

		print('	Ok:', con)
		print('	kO:', nocon)
		print()

	return netcons
	

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
	list_spike = []
	for cell in cells:
		list_v.append(neuron.h.Vector().record(cell.soma[0](0.5)._ref_v))

		spike_times = neuron.h.Vector()
		_spike_detector = neuron.h.NetCon(cell.soma[0](0.5)._ref_v, None, sec=cell.soma[0])
		_spike_detector.record(spike_times)
		list_spike.append(spike_times)

	time_v = neuron.h.Vector().record(neuron.h._ref_t)
	run_run_run(config['tstop'])

	d_result = {'t':list(time_v), 'template_cell_ids':template_cell_ids, 'record_soma_v':{}, 'record_spike_times':{}}

	for num_cell in range(len(list_v)):
		d_result['record_soma_v'][template_cell_ids[num_cell]] = list(list_v[num_cell])
		d_result['record_spike_times'][template_cell_ids[num_cell]] = list(list_spike[num_cell])
	
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

	config = {'tstop':150, 'iclamp_delay':50, 'iclamp_dur':50, 'iclamp_amp':0.1}
	template_cell_ids, R_syn_pre_cell_type, syn_pre_cell_type, num_pre_syn_type, dict_with_info = get_lists(save_dict=False)
	print(template_cell_ids)
	cells, cells_map = load_cells(template_cell_ids)

	# plotting(cells)
	# check_axon_voltage(cells)
	# change_ion_conc(cells, syn_pre_cell_type)
	synapses_map = _create_synapses(gid_cid=template_cell_ids,
									gid_num_pre_syn_mtype=num_pre_syn_type, 
									gid_syn_pre_cell_type=syn_pre_cell_type, 
									gid_R_syn_pre_cell_type=R_syn_pre_cell_type)

	netcons = connect_neurons(synapses_map, cells_map, template_cell_ids, syn_pre_cell_type, False)
	set_stimuls_and_go(config, cells, template_cell_ids, netcons)
	