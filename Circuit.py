import os
import numpy
import json
from neuron import h

### MPI must be initialized before we create a ParallelContext object
h.nrnmpi_init()
pc = h.ParallelContext()

h.load_file("stdrun.hoc")
h.load_file("import3d.hoc")

class Cell(object):
	def __init__(self, gid, cell_id, pc_id, gid_list):
		self._cell_id = cell_id
		self._gid = gid
		self._pc_id = pc_id
		self._gid_list = gid_list
		self._ncs = []

		h.load_file(os.path.join('template', f'{cell_id}.hoc'))
		self._cell = getattr(h, cell_id)(1)
		self._soma = self._cell.soma[0]
		self._synapses = self._cell.synapses
		self._all = self._soma.wholetree()

		self._cell.soma_v = h.Vector().record(self._soma(0.5)._ref_v)
		self._spike_detector = h.NetCon(self._soma(0.5)._ref_v, None, sec=self._soma)
		self._spike_times = h.Vector()
		self._spike_detector.record(self._spike_times)

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

		# if pc.gid_exists(3):
		# 	self.iclamp = h.IClamp(pc.gid2cell(3).soma[0](0.5))
		# 	self.iclamp.delay = 100
		# 	self.iclamp.dur = 300
		# 	self.iclamp.amp = 0.1

		self.l_iclamp = []
		for gid in self.gidlist:
			if pc.gid_exists(gid):
				print(f'		gid_exists{gid}')
				iclamp = h.IClamp(pc.gid2cell(gid).soma[0](0.5))
				iclamp.delay = self.config['iclamp_delay']
				iclamp.dur = self.config['iclamp_dur']
				iclamp.amp = self.config['iclamp_amp']
				self.l_iclamp.append(iclamp)

	def _get_synapses_pre_mtypes(self, filename):
		d_mtype_map, d_synapses, d_type_synapses = {}, {}, {}

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
					d_synapses[synapse_id] = str_mtype
					d_type_synapses[synapse_id] = excitatory
					d_num_mtypes_synapses.setdefault(str_mtype, 0)
					d_num_mtypes_synapses[str_mtype] += 1

		return d_num_mtypes_synapses, d_synapses, d_type_synapses


	def _get_lists(self, save_synapses):	
		self.dict_with_info = {}
		self.syn_pre_cell_type, self.template_cell_ids = [], []
		for _template_paths in os.listdir('template'):
			if _template_paths.endswith('.hoc'):
				cell_id = _template_paths.split('.hoc')[0]
				
				
				num_mtypes_synapses, mtype_synapses, ie_type_synapses = self._get_synapses_pre_mtypes(os.path.join('synapses', cell_id))
				if save_synapses:
					self.dict_with_info[cell_id] = {'num_synapses':len(mtype_synapses), 
					'num_mtypes_synapse':num_mtypes_synapse,
					'mtype_synapses':mtype_synapses, 
					'ie_type_synapses':ie_type_synapses}

				self.syn_pre_cell_type.append(mtype_synapses)
				self.template_cell_ids.append(cell_id)
		if save_synapses:
			with open('synapse_info.json', 'w') as outfile:
				json.dump(self.dict_with_info, outfile)

		return self.template_cell_ids, self.syn_pre_cell_type, self.dict_with_info

	def _set_gids(self):
		#set gids per host
		#each host has an unique id pc.id() from 0 to pc.nhost() - 1.
		self.gidlist = list(range(pc.id(), self.num_neurons, pc.nhost()))
		print(f'\npc.id:{pc.id()} || gidlist:{self.gidlist}, N:{self.num_neurons}, pc.nhost:{pc.nhost()}')
		for gid in self.gidlist:
			# print(f'	gid:{gid}, pc.id:{pc.id()}')
			pc.set_gid2node(gid, pc.id())
	
	def _load_cells(self):
		self.cells = []
		left_load = len(self.gidlist)
		for gid in self.gidlist:#only create the cells that exist on this host
			cell_id = self.template_cell_ids[gid]
			self.cells.append(Cell(gid, cell_id, pc.id(), self.gidlist))

			print(f'	pc.id:{pc.id()} gid:{gid}  cell_id:{cell_id} | left_load:{left_load}')
			left_load -= 1

		### associate the cell with this host and gid
		for cell in self.cells:
			# print(cell.gid, cell.cell_id, cell.soma_v)
			pc.cell(cell._gid, cell._spike_detector)

	def _get_short_inf(self, inf, s_list):
		res = {'L1':0, 'L23':0,'L4':0,'L5':0,'L6':0}
		for syn in s_list:
			layer = inf[int(syn.synapseID)].split('_')[0]
			res[layer] += 1
		return res

	def _connect_cells(self):
		self.log_pre_post_syn = {}
		for gid_t, target in zip(self.gidlist, self.cells):
			
			post_syn_cell_id = target._cell_id
			self.log_pre_post_syn[post_syn_cell_id] = {}
			info_for_target = self.syn_pre_cell_type[gid_t]#{target_syn_id: 'L1_HAC', 1: 'L1_HAC', 2:'L2_DSF'}
			_inf = self._get_short_inf(info_for_target, target._synapses.synapse_list)
			self.log_pre_post_syn[post_syn_cell_id]['info_for_target'] = _inf
			print('post_syn:', target, '	num pre_syn:', len(target._synapses.synapse_list), _inf)
			
			for gid_s in range(self.num_neurons):
				if gid_t == gid_s:
					continue

				pre_syn_cell_id = self.template_cell_ids[gid_s]
				pre_syn_layer = pre_syn_cell_id.split('_')[0]

				if not pre_syn_cell_id in self.log_pre_post_syn[post_syn_cell_id]:
					self.log_pre_post_syn[post_syn_cell_id][pre_syn_cell_id] = 0
				# print('	pre_syn:', pre_syn_cell_id, end = '->	')
				_cnt=0
				for syn in target._synapses.synapse_list:
					could_be_pre_syn = info_for_target[int(syn.synapseID)]
					could_be_pre_syn_layer = could_be_pre_syn.split('_')[0]

					if could_be_pre_syn_layer == pre_syn_layer:
						# print(could_be_pre_syn, end=' ')
						# if _cnt > 10:
						# 	_cnt = 0
						# 	print('\n			', end='')

						nc = pc.gid_connect(gid_s, syn)
						nc.weight[0] = self.syn_w
						nc.delay = self.syn_delay
						target._ncs.append(nc)
						self.log_pre_post_syn[post_syn_cell_id][pre_syn_cell_id] += 1
						_cnt+=1
					# print()
			print('all presyn connected to postsyn:', target, end='\n\n')

