import matplotlib as mpl
import matplotlib.pyplot as plt
import json
import os
import numpy as np
import time
import random
import seaborn as sns

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
	plt.savefig('Fig30.png')


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
			return i
	print(cellid)
	exit('get_mtype')

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
	if not os.path.exists(opj(res_folder, 'Mpre_syn.json')):
	# if os.path.exists(opj(res_folder, 'Mpre_syn.json')):
		dict_connections = {}
		# file = 'synapses_0_map.json'	
		# file = 'synapses_1_map.json'
		file = 'synapses_3_map.json'
		
		with open(opj(res_folder, file), 'r') as outfile:
			result = json.load(outfile)
		Mpre_syn={}
		Mpost_syn={}
		for k in excit.keys():
			Mpre_syn[k]=[]#how many uniq preneurons cell id in this mtype
			Mpost_syn[k]=[]

		all_m_pre = []
		lll = len(result.keys())
		for post_syn, pre_synaptics in result.items():
			print(lll)
			lll -= 1
			mtype_post_syn = get_mtype(post_syn)
			
			uniq = list(np.unique(list(pre_synaptics.values())))

			Mpre_syn[mtype_post_syn].extend(uniq)

			dpre_syn = {pre_syn:0 for pre_syn in uniq}
			
			for num, pre_syn in pre_synaptics.items():
				dpre_syn[pre_syn] += 1

			for pre_syn in uniq:
				mtype_pre_syn = get_mtype(pre_syn)
				if inhib[mtype_pre_syn] == 1:
					im[mtype_post_syn] += dpre_syn[pre_syn]
				else:
					em[mtype_post_syn] += dpre_syn[pre_syn]

				Mpost_syn[mtype_pre_syn].append(post_syn)

		for k in excit.keys():
			print(k)
			Mpost_syn[k] = len(np.unique(Mpost_syn[k]))
			Mpre_syn[k] = len(np.unique(Mpre_syn[k]))

		with open(opj(res_folder, 'em.json'), 'w') as outfile:
			json.dump(em, outfile)
		with open(opj(res_folder, 'im.json'), 'w') as outfile:
			json.dump(im, outfile)

		with open(opj(res_folder, 'Mpre_syn.json'), 'w') as outfile:
			json.dump(Mpre_syn, outfile)
		with open(opj(res_folder, 'Mpost_syn.json'), 'w') as outfile:
			json.dump(Mpost_syn, outfile)
	else:
		with open(opj(res_folder, 'em.json'), 'r') as outfile:
			em = json.load(outfile)
		with open(opj(res_folder, 'im.json'), 'r') as outfile:
			im = json.load(outfile)

		# res_folder = '../Connectome'
		with open(opj(res_folder, 'Mpost_syn.json'), 'r') as outfile:
			Mpost_syn = json.load(outfile)
		with open(opj(res_folder, 'Mpre_syn.json'), 'r') as outfile:
			Mpre_syn = json.load(outfile)
	
	
	
	# mpl.rcParams['font.size'] = 33.0
	with open('BBPjson/pathways_anatomy_factsheets_simplified.json', 'r') as outfile:
		als = json.load(outfile)

	#################
	# for prepost, v in als.items():
	# 	vvv = v['total_synapse_count']
	# 	prepost = prepost.replace('-', '_')
	# 	pre, post = prepost.split(':')
	# 	if inhib[pre] == 1:
	# 		im[post] += v['total_synapse_count']
	# 	else:
	# 		em[post] += v['total_synapse_count']
	#################

	# # plt.figure(figsize=(10, 10))
	# ipost, ipre = {}, {}
	# epost, epre = {}, {}
	# for k, v in excit.items():
	# 	if v != 0:
	# 		epost[k], epre[k] = 0,0
	# for k, v in inhib.items():
	# 	if v != 0:
	# 		ipost[k], ipre[k] = 0,0
	# for prepost, v in als.items():
	# 	# vvv = v['number_of_divergent_neuron_std'] + v['number_of_convergent_neuron_mean']
	# 	vvv = v['total_synapse_count']
	# 	prepost = prepost.replace('-', '_')
	# 	pre, post = prepost.split(':')
	# 	if inhib[pre] == 1:
	# 		ipost[pre]+=vvv
	# 	else:
	# 		epost[pre]+=vvv

	# 	if inhib[post] == 1:
	# 		ipre[post]+=vvv
	# 	else:
	# 		epre[post]+=vvv


	Mpre_syn = {k.replace('-','_'):v for k, v in Mpre_syn.items()}
	Mpost_syn = {k.replace('-','_'):v for k, v in Mpost_syn.items()}
	ipost, ipre = {}, {}
	epost, epre = {}, {}
	for k, v in excit.items():
		if v != 0:
			epost[k], epre[k] = Mpost_syn[k], Mpre_syn[k]
			print(epost[k], epre[k], '[k]')
	for k, v in inhib.items():
		if v != 0:
			ipost[k], ipre[k] = Mpost_syn[k], Mpre_syn[k]
			print(ipost[k], ipre[k], 'ipost[k], ipre[k]')

	epre = get_reverse(epre)
	epost = get_reverse(epost)
	ipre = get_reverse(ipre)
	ipost = get_reverse(ipost)

	# ax = plt.subplot(2,2,1, frameon=False)
	# ax.barh(list(epost.keys()), list(epre.values()), color='red')#, height=0.4, color='blue')
	# ax.set_xlim([max(list(epre.values())), 0])  # decreasing time
	# ax.set_yticklabels([])
	# ax = plt.subplot(2,2,2, frameon=False)
	# ax.barh(list(epost.keys()), list(epost.values()), color='blue')#, height=0.4, color='blue')

	# # ax = plt.subplot(1,2,1, frameon=False)
	# ax = plt.subplot(2,2,3, frameon=False)
	# ax.barh(list(ipost.keys()), list(ipre.values()), color='red')#, height=0.4, color='blue')
	# ax.set_xlim([max(list(ipre.values())), 0])  # decreasing time
	# ax.set_yticklabels([])

	# # ax = plt.subplot(1,2,2, frameon=False)
	# ax = plt.subplot(2,2,4, frameon=False)
	# ax.barh(list(ipost.keys()), list(ipost.values()), color='blue')#, height=0.4, color='blue')
	# print(sum(list(ipost.values()))+ sum(list(epost.values())), sum(list(ipre.values()))+ sum(list(epre.values())))
	# plt.show()

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
	plt.savefig('fig32.png')

	


def pplot_per_layer(template_cell_ids, t, result, t1=0,t2=0, folder):

	template_cell_ids.reverse()
	if not (template_cell_ids[0].startswith('L6') and template_cell_ids[-1].startswith('L1')):
		print(template_cell_ids[0], template_cell_ids[-1])
		exit('-')

	mpl.rcParams['font.size'] = 15.0
	l = len(template_cell_ids)
	print('pplot_per_layer', len(template_cell_ids))
	print(template_cell_ids[0], template_cell_ids[-1])
	plt.figure(figsize=(14,14))
	
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


def pvlines_per_layer(template_cell_ids, t, result, t1=0,t2=0, figsize=(13,9), folder):
	global should_mtypes
	
	template_cell_ids.reverse()
	if not (template_cell_ids[0].startswith('L6') and template_cell_ids[-1].startswith('L1')):
		print(template_cell_ids[0], template_cell_ids[-1])
		exit('-')
	mpl.rcParams['font.size'] = 15.0
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

def hihistogramm(t1,t2,t, all_spike_histogramm, folder):
	if t1 !=0 and t2 != 0:
		b = np.linspace(t1,t2, 100, endpoint=True)
	else:
		b = np.linspace(0,max(t), 500, endpoint=True)


	spike_histogramm = {'L1':[], 'L23':[], 'L4':[], 'L5':[], 'L6':[]}
	for k,v in all_spike_histogramm.items():
		spike_histogramm[k.split('_')[0]].extend(v)
	l = len(spike_histogramm)
	plt.figure(figsize=(16,9))
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
	plt.savefig(opj(folder, '1hihistogramm.png'))

def hihistogramm3(t1,t2,t, spike_histogramm, folder):
	if t1 !=0 and t2 != 0:
		b = np.linspace(t1,t2, 100, endpoint=True)
	else:
		b = np.linspace(0,max(t), 250, endpoint=True)

	l = len(spike_histogramm)
	plt.figure(figsize=(16,9))
	for i, (cell_id, m_spikes) in enumerate(spike_histogramm.items()):
		ax = plt.subplot(l, 1, i+1, frameon=False)
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
	# plt.show()
	plt.savefig(opj(folder, '1hihistogramm3.png'))

def plot_parallel_cores(t1, t2, res_folder='result'):
	global should_mtypes, should_etypes
	# inhib, excit = get_inh_exi()
	# create_fig30(inhib, excit, res_folder)
	# exit()
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
				frstt = -1
				lstt = -1
				if t1 > 0 and t2 > 0:	
					for i in range(len(t)):
						if frstt == -1 and t1 < t[i]:
							frstt = i
						if lstt == -1 and t[i] > t2:
							lstt = i
							break
					t = [i for i in t if t1 < i < t2]

			for iii, (cell_id, soma_v) in enumerate(result['record_soma_v'].items()):
				spiket = result['record_spike_times'][cell_id]
				for m_cell_id in should_mtypes:
					if cell_id.startswith(m_cell_id):
						break
				
				layer_nums[cell_id.split('_')[0]] += 1
				if t1 > 0 and t2 > 0:
					spiket = [i for i in spiket if t1 < i < t2]

				spike_times.setdefault(m_cell_id, [])
				spike_times[m_cell_id].extend(spiket)

				cell_spike_times[cell_id] = np.asarray(spiket)

				isi = []
				if len(spiket) > 1:
					# isi = [r-l for r,l in zip(spiket[1:], spiket[:-1])]
					isi = spiket
					hisst[m_cell_id].extend(isi)
				
				soma_v = result['record_soma_v'][cell_id]
				if frstt != -1 and lstt != -1:
					soma_v = soma_v[frstt:lstt]
				soma_v = np.asarray(soma_v)
				
				if not m_cell_id in short_soma_v.keys():
					short_soma_v[m_cell_id] = soma_v
				# else:
				# 	was_soma_v = short_soma_v[m_cell_id]
				# 	if np.std(soma_v) > np.std(was_soma_v):
				# 		short_soma_v[m_cell_id] = soma_v
				
				num_neurons += 1

				# if len(cell_spike_times.keys()) != num_neurons:
				# 	print(cell_id)
				# 	exit('---1')
			# if len(short_soma_v.keys()) > 50:
			# 	break
			# break
			# if num_neurons > 400:
			# 	break

	print(layer_nums)
	print(sum(layer_nums.values()), 'SUM')
	print(len(cell_spike_times.keys()), num_neurons, '----------------')

	for cell_id in spike_times.keys():
		spike_times[cell_id] = np.asarray(spike_times[cell_id])

	step = 10
	for window in range(t1-10, t2, step):
		layer_activity = {'L1':0, 'L23':0, 'L4':0, 'L5':0, 'L6':0}
		for cell_id, spiket in spike_times.items():
			layer = cell_id.split('_')[0]
			layer_activity[layer] += len(spiket[(spiket > window) & (spiket < window + 10)])
		if sum(layer_activity.values())>0:
			layer_activity = dict(sorted(layer_activity.items(), key=lambda x: x[1], reverse=True))
			print(f'{window}:{window + step} 	{layer_activity}')

	sorted_cell_spike, wascellid = get_sorted_cellids(list(cell_spike_times.keys()), should_mtypes)
	sorted_short_soma, wassoma = get_sorted_cellids(list(short_soma_v.keys()), should_mtypes)
	if not os.path.exists(opj(res_folder, 'cell.json')):
		result = {'cell_ids':sorted_cell_spike}
		with open(opj(res_folder, 'cell.json'), 'w') as fd:
			json.dump(result, fd)

	print('nums:', wascellid, wascellid==len(sorted_cell_spike)==len(np.unique(sorted_cell_spike)),  'num_neurons:', num_neurons)
	hihistogramm3(t1,t2,t, hisst)
	hihistogramm(t1,t2,t, hisst)
	figsize=(7,14)
	if t1 == 0 and t2 == 0:
		figsize=(14,14)
	pvlines_per_layer(sorted_cell_spike, t, cell_spike_times, t1,t2, figsize=figsize)
	pplot_per_layer(sorted_short_soma, t, short_soma_v, t1, t2)


	


if __name__ == '__main__':
	should_mtypes = ['L1_DAC','L1_NGC_DA','L1_NGC_SA','L1_HAC','L1_DLAC','L1_SLAC','L23_PC','L23_MC','L23_BTC','L23_DBC','L23_BP','L23_NGC','L23_LBC','L23_NBC','L23_SBC','L23_ChC','L4_PC','L4_SP','L4_SS','L4_MC','L4_BTC','L4_DBC','L4_BP','L4_NGC','L4_LBC','L4_NBC','L4_SBC','L4_ChC','L5_TTPC1','L5_TTPC2','L5_UTPC','L5_STPC','L5_MC','L5_BTC','L5_DBC','L5_BP','L5_NGC','L5_LBC','L5_NBC','L5_SBC','L5_ChC','L6_TPC_L1','L6_TPC_L4','L6_UTPC','L6_IPC','L6_BPC','L6_MC','L6_BTC','L6_DBC','L6_BP','L6_NGC','L6_LBC','L6_NBC','L6_SBC','L6_ChC']
	should_etypes = ['cAC', 'bAC', 'cNAC', 'bNAC', 'dNAC', 'cSTUT', 'bSTUT', 'dSTUT', 'cIR' , 'bIR', 'cADpyr']

	opj = os.path.join
	t1, t2 = 130, 370
	# t1, t2 = 0, 0
	

	# plot_parallel_cores(t1,t2, res_folder='../CloudResult2.5sec_100percent_2mmol')
	# plot_synapses(res_folder='../CorrectCloudResult2.5sec_100percent_2mmol')
	# exit()

	# plot_parallel_cores(t1,t2, res_folder='result')
	# plot_parallel_cores(t1,t2, res_folder='../CloudResult2.5sec_95percent_2mmol'
	plot_parallel_cores(t1,t2,res_folder='../CorrectCloudResult2.5sec_100percent_2mmol')
	
	