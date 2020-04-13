import matplotlib.pyplot as plt
import json
import os
import numpy as np

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


def pplot_per_layer(template_cell_ids, t, result):
	l = len(template_cell_ids)
	plt.figure(figsize=(9,9))
	for i, cell_id in enumerate(template_cell_ids):
		ax = plt.subplot(l, 1, i+1, frameon=False)
		plt.plot(t, result[cell_id])
		plt.ylim([-90,50])
		h = ax.set_ylabel(cell_id)
		h.set_rotation(0)
		if i + 1 == l:
			h = ax.set_xlabel('t')
		else:
			ax.set_xticklabels([])
			ax.set_yticklabels([])

	plt.show()


def pvlines_per_layer(template_cell_ids, t, result):
	l = len(template_cell_ids)
	plt.figure(figsize=(9,9))
	for i, cell_id in enumerate(template_cell_ids):
		ax = plt.subplot(l, 1, i+1, frameon=False)
		# ax = plt.subplot(1, 1, 1)
		for spike in result[cell_id]:
			plt.vlines(spike, 0, 1)
		plt.ylim([0,1])
		plt.xlim([0,max(t)])
		h = ax.set_ylabel(cell_id)
		h.set_rotation(0)
		if i + 1 == l:
			h = ax.set_xlabel('t')
		else:
			ax.set_xticklabels([])
			ax.set_yticklabels([])

	plt.show()


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

		
def plot_parallel_cores():
	template_cell_ids, t = [], []
	res_soma_v = {}
	res_spike_times = {}

	for file in os.listdir('result'):
		if file.startswith('host_') and file.endswith('.json'):
			print(file)
			with open(f'result/{file}', 'r') as outfile:
				result = json.load(outfile)
			if not len(t):
				t = result['t']
			
			for cell_id, soma_v in result['record_soma_v'].items():
				res_soma_v[cell_id] = soma_v
				res_spike_times[cell_id] = result['record_spike_times'][cell_id]
				print('	', cell_id, max(soma_v))
				template_cell_ids.append(cell_id)
	template_cell_ids.sort()
	pvlines_per_layer(template_cell_ids, t, res_spike_times)
	exit()
	if len(template_cell_ids) > 100:
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
	
	pplot_per_layer(template_cell_ids, t, res_soma_v)

if __name__ == '__main__':
	# plot_single_core()
	# plot_all_cells_parallel()
	plot_parallel_cores()
