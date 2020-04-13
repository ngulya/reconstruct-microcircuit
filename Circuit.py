import os
import numpy
import json
from neuron import h
import random

random.seed(13)

### MPI must be initialized before we create a ParallelContext object
h.nrnmpi_init()
pc = h.ParallelContext()

h.load_file("stdrun.hoc")
h.load_file("import3d.hoc")

class Cell(object):
	def __init__(self, gid, cell_id, pc_id, gid_host_list):
		self._cell_id = cell_id#'L1_DAC_bSyn_192_5'#unique
		self._gid = gid#integer unique
		self._pc_id = pc_id#host id
		self._gid_host_list = gid_host_list
		self._ncs = []

		h.load_file(os.path.join('template', f'{cell_id}.hoc'))
		self._cell = getattr(h, cell_id)(1)
		self._soma = self._cell.soma[0]
		self._synapses = self._cell.synapses
		self._all = self._soma.wholetree()

		self.soma_v = h.Vector().record(self._soma(0.5)._ref_v)
		self.spike_times = h.Vector()

		self._spike_detector = h.NetCon(self._soma(0.5)._ref_v, None, sec=self._soma)
		self._spike_detector.record(self.spike_times)

	def __repr__(self):
		return f'Cell[{self._gid}]_{self._cell}'

class Circuit(object):

	def __init__(self, config, num_neurons=5, syn_w=0.05, syn_delay=1, save_synapses=False):
		self.config = config
		self.num_neurons = num_neurons
		self.syn_w = syn_w
		self.syn_delay = syn_delay
		self._get_lists(save_synapses)
		self._set_gids()### assign gids to processors
		self._load_cells()
		self._connect_cells()
		self._set_stimuls()

	def _get_synapses_pre_mtypes(self, filename):
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


	def _get_lists(self, save_synapses):	
		self.dict_with_info = {}
		self.cid_gid = {}

		self.gid_num_pre_syn_mtype, self.gid_cid = [], []
		self.gid_syn_pre_cell_type, self.gid_R_syn_pre_cell_type = [], []
		for _template_paths in os.listdir('template'):
			if _template_paths.endswith('.hoc'):
				cell_id = _template_paths.split('.hoc')[0]
				
				
				num_mtypes_synapses, R_mtype_synapses, mtype_synapses, ie_type_synapses = self._get_synapses_pre_mtypes(os.path.join('synapses', cell_id))
				'''
				synapse_id = [0...100...N]

				num_mtypes_synapses[m-type] = 2 		|nums pre_synaptic cell
				mtype_synapses[0..N] = 'L1_DAC'.. 		|m-type
				R_mtype_synapses['L1_DAC'] = [0,1,10]	|m-type
				ie_type_synapses[0..N] = 0/1  			|0-inhibitory, 1-excitatory
				'''
				if save_synapses:
					self.dict_with_info[cell_id] = {'num_synapses':len(mtype_synapses), 
					'num_mtypes_synapses':num_mtypes_synapses,
					'mtype_synapses':mtype_synapses, 
					'ie_type_synapses':ie_type_synapses}

				
				self.cid_gid[cell_id] = len(self.gid_cid)
				self.gid_cid.append(cell_id)
				self.gid_num_pre_syn_mtype.append(num_mtypes_synapses)
				self.gid_syn_pre_cell_type.append(mtype_synapses)
				self.gid_R_syn_pre_cell_type.append(R_mtype_synapses)

		if save_synapses:
			with open('synapse_info.json', 'w') as outfile:
				json.dump(self.dict_with_info, outfile)

	def _set_gids(self):
		#set gids per host
		#each host has an unique id pc.id() from 0 to pc.nhost() - 1.
		self.gid_host_list = list(range(pc.id(), self.num_neurons, pc.nhost()))
		print(f'\npc.id:{pc.id()} || gidlist:{self.gid_host_list}, N:{self.num_neurons}, pc.nhost:{pc.nhost()}')
		for gid in self.gid_host_list:
			# print(f'	gid:{gid}, pc.id:{pc.id()}')
			pc.set_gid2node(gid, pc.id())
	
	def _load_cells(self):
		self.obj_cells = []
		self.obj_cell_cid = {}
		left_load = len(self.gid_host_list)
		for gid in self.gid_host_list:#only create the cells that exist on this host
			cell_id = self.gid_cid[gid]
			_cell = Cell(gid, cell_id, pc.id(), self.gid_host_list)

			self.obj_cells.append(_cell)
			self.obj_cell_cid[cell_id] = _cell

			print(f'	pc.id:{pc.id()} gid:{gid}  cell_id:{cell_id} | left_load:{left_load}')
			left_load -= 1

		### associate the cell with this host and gid
		for cell in self.obj_cells:
			# print(cell.gid, cell.cell_id, cell.soma_v)
			pc.cell(cell._gid, cell._spike_detector)

	def _sort_dict(self, d):
		res_d = {}
		list_key = list(d.keys())
		list_key.sort()
		for k in list_key:
			res_d[k] = d[k]
		return res_d

	def _get_logs(self, gid_t, s_list):
		short_log, full_log = {}, {}

		for syn in s_list:
			pre_syn_mtype = self.gid_syn_pre_cell_type[gid_t][int(syn.synapseID)]#{synapseID1: 'L1_HAC', synapseID2: 'L1_HAC'}
			layer = pre_syn_mtype.split('_')[0]

			short_log.setdefault(layer, 0)
			short_log[layer] += 1

			full_log.setdefault(pre_syn_mtype, 0)
			full_log[pre_syn_mtype] += 1

		return self._sort_dict(short_log), self._sort_dict(full_log)

	def _has_this_mtype(self, target_num, mtype):
		result = []
		for _target_num in range(len(self.gid_cid)): 
			if mtype in self.gid_cid[_target_num] and _target_num != target_num:
				result.append(_target_num)
				
		return result

	def _create_synapses(self):
		# self.gid_cid[0]= 'L1_HAC_bNAC219_3'
		# self.gid_num_pre_syn_mtype[0]= {'L1_HAC': 19, 'L23_NBC': 32, 'L23_PC': 21, 'L23_MC': 69,'L6_NBC_TP': 3}
		# self.gid_syn_pre_cell_type[0]{0:'L1_DAC'}.. 		|m-type
		# self.gid_R_syn_pre_cell_type[0]{'L1_DAC': [0,1,3,4]}
		synapses_map = {}
		for target_num in range(len(self.gid_cid)):
			target_cell_id = self.gid_cid[target_num]
			pre_syn_mtype = self.gid_num_pre_syn_mtype[target_num]
			mtype_synapses = self.gid_syn_pre_cell_type[target_num]
			R_mtype_synapses = self.gid_R_syn_pre_cell_type[target_num]
			
			synapses_map[target_cell_id] = mtype_synapses.copy()#{0:'L1_DAE'}

			for mtype, num_synapses in pre_syn_mtype.items():
				available = self._has_this_mtype(target_num, mtype)
				num_available = len(available)
				
				if num_available:
					drop = False
					if num_synapses <= num_available:
						drop = True
					
					for syn_id in R_mtype_synapses[mtype]:
						cell_id = random.choice(available)
						if drop:
							available.remove(cell_id)
						synapses_map[target_cell_id][syn_id] = self.gid_cid[cell_id]
				else:
					for syn_id in R_mtype_synapses[mtype]:
						synapses_map[target_cell_id][syn_id] = 'no-data'

		return synapses_map

	def _connect_cells(self, _create_synapses=True, synapse_json=None):
		# print('self.gid_cid[]=', self.gid_cid)
		# print()
		# print('self.gid_num_pre_syn_mtype[0]=', self.gid_num_pre_syn_mtype[0])
		synapses_map = self._create_synapses()
		# exit('---')
		self.log_pre_post_syn = {}
		for gid_t, target in zip(self.gid_host_list, self.obj_cells):
			
			post_syn_cid = target._cell_id
			self.log_pre_post_syn[post_syn_cid] = {}
			
			_short_log, _full_log = self._get_logs(gid_t, target._synapses.synapse_list)
			self.log_pre_post_syn[post_syn_cid]['info_for_target'] = _short_log
			self.log_pre_post_syn[post_syn_cid]['info_for_target_full'] = _full_log

			print('post_syn:', target, '	num pre_syn:', len(target._synapses.synapse_list), _short_log)
			
			for syn_id, pre_syn_cid in synapses_map[post_syn_cid].items():
				if pre_syn_cid == 'no-data':
					continue
				synapse = target._synapses.synapse_list[syn_id]
				could_be_pre_syn = self.gid_syn_pre_cell_type[gid_t][int(synapse.synapseID)]
				if not could_be_pre_syn in pre_syn_cid:
					exit(f'Synapses Error {could_be_pre_syn} != {pre_syn_cid}')

				pre_syn_gid = self.cid_gid[pre_syn_cid]

				nc = pc.gid_connect(pre_syn_gid, synapse)
				nc.weight[0] = self.syn_w
				nc.delay = self.syn_delay
				target._ncs.append(nc)

				self.log_pre_post_syn[post_syn_cid].setdefault(pre_syn_cid, 0)
				self.log_pre_post_syn[post_syn_cid][pre_syn_cid] += 1

		return 1

	def _set_stimuls(self):
		# if pc.gid_exists(3):
		# 	self.iclamp = h.IClamp(pc.gid2cell(3).soma[0](0.5))
		# 	self.iclamp.delay = 100
		# 	self.iclamp.dur = 300
		# 	self.iclamp.amp = 0.1
		self.l_iclamp = []
		for gid in self.gid_host_list:
			if pc.gid_exists(gid):
				print(f' On host {pc.id()}	gid {gid} exists as {self.gid_cid[gid]}')
				iclamp = h.IClamp(pc.gid2cell(gid).soma[0](0.5))
				iclamp.delay = self.config['iclamp_delay']
				iclamp.dur = self.config['iclamp_dur']
				iclamp.amp = self.config['iclamp_amp']
				self.l_iclamp.append(iclamp)


