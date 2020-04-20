import matplotlib.pyplot as plt
import json
import os
import numpy as np
import time
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


def pplot_per_layer(template_cell_ids, t, result, t1=0,t2=0):
	l = len(template_cell_ids)
	plt.figure(figsize=(9,9))
	for i, cell_id in enumerate(template_cell_ids):
		ax = plt.subplot(l, 1, i+1, frameon=False)
		plt.plot(t, result[cell_id])
		if t1 != 0 and t2 != 0:
			plt.axvline(t1, -90, 50, color='green')
			plt.axvline(t2, -90, 50, color='green')

		plt.ylim([-90,50])
		h = ax.set_ylabel(cell_id)
		h.set_rotation(0)
		if i + 1 == l:
			h = ax.set_xlabel('t')
		else:
			ax.set_xticklabels([])
			ax.set_yticklabels([])

	plt.show()


def pvlines_per_layer(template_cell_ids, t, result, t1=0,t2=0):
	

	# l = len(template_cell_ids)
	# template_cell_ids.reverse()
	# cell_names = {'L1':0, 'L23':0, 'L4':0, 'L5':0, 'L6':0}
	# for i in template_cell_ids:
	# 	cell_names[i.split('_')[0]] += 1
	# fig, ax = plt.subplots(figsize=(13,9))
	# for ti in range(0, int(float(max(t))/10)*10, 10):
	# 	plt.vlines(ti, 0, l+1, color='blue', linewidths=1)
	# for i, cell_id in enumerate(template_cell_ids):
	# 	for spike in result[cell_id]:
	# 		plt.vlines(spike, i, i+1)
	# if t1 != 0 and t2 != 0:
	# 	plt.vlines(t1, 0, l+1, color='green')
	# 	plt.vlines(t2, 0, l+1, color='green')
	# ax.set_yticklabels(['L6', 'L5', 'L4', 'L23', 'L1'])
	# ax.set_xlabel(cell_names)
	# plt.ylim([0,l+1])
	# plt.xlim([0,max(t)])
	# plt.show()
	l = len(template_cell_ids)
	plt.figure(figsize=(16,9))
	for i, cell_id in enumerate(template_cell_ids):
		ax = plt.subplot(l, 1, i+1, frameon=False)
		# ax = plt.subplot(1, 1, 1)
		for spike in result[cell_id]:
			plt.vlines(spike, 0, 1)
		
		if t1 != 0 and t2 != 0:
			plt.axvline(t1, -90, 50, color='green')
			plt.axvline(t2, -90, 50, color='green')

		plt.ylim([0,1])
		plt.xlim([0,max(t)])
		if i%5 == 0:
			cell_id = '_'.join(cell_id.split('_')[:3])
			h = ax.set_ylabel(cell_id)
			h.set_rotation(0)
		if i < l -1:
			ax.set_xticklabels([])
			ax.set_yticklabels([])

	h = ax.set_xlabel('t')
	plt.show()
	exit()


def plot_single_core():
	if os.path.exists('result/single.json'):
		with open('result/single.json', 'r') as outfile:
			result = json.load(outfile)

		template_cell_ids = result['template_cell_ids']
		template_cell_ids.sort()
		#result['record_soma_v'] = {'L1_Dac_121_4':[],.....}
		#template_cell_ids = ['L1...', 'L2....']
		pplot_per_layer(template_cell_ids, result['t'], result['record_soma_v'])
		pvlines_per_layer(template_cell_ids, result['t'], result['record_spike_times'])

		
def plot_parallel_cores(t1,t2, tanalyze):
	template_cell_ids, t = [], []
	res_soma_v = {}
	res_spike_times = {}
	num_neurons = 0
	for file in os.listdir('result'):
		if file.startswith('host_') and file.endswith('.json'):
			print(file)
			with open(f'result/{file}', 'r') as outfile:
				result = json.load(outfile)
			if not len(t):
				t = result['t']
			
			
			most_spiked = {}
			for cell_id, spike_times in result['record_spike_times'].items():
				# 'L6_IPC_cADpyr231_3_dplc_68'
				# print(cell_id)
				m_cell_id = cell_id
				# m_cell_id = '_'.join(cell_id.split('_')[:-4])
				# print(m_cell_id)
				# exit()
				# 'L6_IPC'
				num_spike = 0
				for sp in spike_times:
					if sp < t1 or sp > t2:
						num_spike += 1
				# most_spiked.setdefault(m_cell_id, [0,''])
				most_spiked.setdefault(m_cell_id, [num_spike,cell_id])
				if most_spiked[m_cell_id][0] < num_spike:
					most_spiked[m_cell_id][0] = num_spike
					most_spiked[m_cell_id][1] = cell_id
				num_neurons += 1
			
			for m_cell_id, list_num_cell_id in most_spiked.items():
				
				cell_id = list_num_cell_id[1]
			
				soma_v = result['record_soma_v'][cell_id]
				spike_times = result['record_spike_times'][cell_id]
				# cell_id = '_'.join(cell_id.split('_')[:-4])
				res_soma_v[cell_id] = np.asarray(soma_v)
				res_spike_times[cell_id] = np.asarray(spike_times)
				
				if list_num_cell_id[0] == 0:
					continue
				print('	', cell_id, max(soma_v))
				template_cell_ids.append(cell_id)

	
	tmax = int(max(t))
	step = 3
	res_spike_times[cell_id] > tanalyze
	for window in range(t2, tmax, step):
		
		layer_activity = {'L1':0, 'L23':0, 'L4':0, 'L5':0, 'L6':0}
		for cell_id, spike_times in res_spike_times.items():
			layer = cell_id.split('_')[0]
			layer_activity[layer] += len(spike_times[(spike_times > window) & (spike_times < window + 10)])
		if sum(layer_activity.values())>0:
			layer_activity = dict(sorted(layer_activity.items(), key=lambda x: x[1], reverse=True))

			print(f'{window}:{window + step} 	{layer_activity}')

	print('nums:', len(template_cell_ids), 'num_neurons:', num_neurons)
	time.sleep(1)
	template_cell_ids.sort()

	if len(template_cell_ids) > 10000:
		print('was:', len(template_cell_ids))
		grouped_mtype = {}
		for cell_id in template_cell_ids:
			mtype = cell_id.split('_')[0] + cell_id.split('_')[1]

			grouped_mtype.setdefault(mtype, [])
			grouped_mtype[mtype].append(res_soma_v[cell_id])
		template_cell_ids = list(grouped_mtype.keys())
		template_cell_ids.sort()
		new_result = {}
		for mtype in template_cell_ids:
			averagel = []
			for i in range(len(t)):
				tmp=[]
				for l in grouped_mtype[mtype]:
					tmp.append(l[i])
				averagel.append(np.mean(tmp))
			new_result[mtype] = averagel
		res_soma_v = new_result
		print('now:', len(template_cell_ids))
	pvlines_per_layer(template_cell_ids, t, res_spike_times, t1,t2)
	pplot_per_layer(template_cell_ids, t, res_soma_v,t1,t2 )

if __name__ == '__main__':
	s = np.asarray([12,34,5])
	t1,t2 = 0,30
	tanalyze = 30
	# plot_single_core()
	# plot_all_cells_parallel()
	plot_parallel_cores(t1,t2, tanalyze)
