from neuron import h
from neuron.units import ms, mV
from Circuit import Circuit
import time
import json
import numpy
import os
# clear; mpiexec -n 6 nrniv 3_run.py

def run_run_run(tstop):
	print(f'\nRunning for {tstop} ms')
	tstart = time.time()
	leftsec = tstop
	num_windows = max([int(float(tstop)/50), 30])
	windows = numpy.linspace(0,tstop,num_windows)

	for iround, neuron.h.tstop in enumerate(windows[1:]):
		print(f'\rLeft {int(tstop-neuron.h.tstop)} ms in simulation || {leftsec} sec', end='    ')
		if iround == 0:
			neuron.h.run()
		else:
			neuron.h.continuerun(neuron.h.tstop)
		tspent = time.time() - tstart
		leftsec = numpy.round(tspent*(tstop-neuron.h.tstop)/neuron.h.tstop, 1)

	Elapsed = numpy.round(time.time() - tstart, 2)
	print(f'\rElapsed time {Elapsed} sec .......')
	return Elapsed

def get_neuron_number(folder='template'):
	num_neurons = 0
	for _template_paths in os.listdir(folder):
		if _template_paths.endswith('.hoc'):
			num_neurons += 1
	return num_neurons

def information_print(CRC):
	log_pre_post_syn = CRC.log_pre_post_syn
	
	for post_syn, d in log_pre_post_syn.items():
		info_for_target_full = log_pre_post_syn[post_syn]['info_for_target_full']
		
		# print(post_syn, info_for_target_full)
		target_num = CRC.cid_gid[post_syn]
		need_pre_syn_mtype = sum(info_for_target_full.values())
		info_for_target_full2 = {}
		for mtype, num_synapses in info_for_target_full.items():
			available = CRC._has_this_mtype(target_num, mtype)
			if len(available):
				info_for_target_full2[mtype] = True
			else:
				info_for_target_full2[mtype] = False

		FLAG=False
		for pre_syn, nums in log_pre_post_syn[post_syn].items():
			if pre_syn in ['info_for_target_full', 'info_for_target']:
				continue
			mtype = 0
			for _mtype, _nums in info_for_target_full.items():
				if _mtype in pre_syn:
					mtype = _mtype
					FLAG = True
			if mtype:
				info_for_target_full[mtype] -= nums
		after_conn_pre_syn_mtype = sum(info_for_target_full.values())

		if after_conn_pre_syn_mtype != 0:
			print(f'\nwas:{need_pre_syn_mtype} now:{after_conn_pre_syn_mtype}!= 0, target:{post_syn} FLAG:{FLAG}')
			for k, v in info_for_target_full.items():
				if v != 0:
					print(f' {k}:{v} ', end = '')
			print()
		for pre_syn, l in d.items():
			if pre_syn in ['info_for_target_full', 'info_for_target']:
				continue
			for name, status in info_for_target_full2.items():
				if name in pre_syn:
					if not status:
						print('	', info_for_target_full2)
						print(f'	{pre_syn}: {l}')
						exit(f'error {pre_syn} | info_for_target_full2')
	
if __name__ == '__main__':
	config = {'tstop':150, 'iclamp_delay':0, 'iclamp_dur':30, 'iclamp_amp':0.1}
	num_neurons = get_neuron_number()

	CRC = Circuit(config, num_neurons)

	pc = h.ParallelContext()
	pc.set_maxstep(10 * ms)

	t = h.Vector().record(h._ref_t)
	h.finitialize(-65 * mV)
	print(f'\nrun №{pc.id()}')# pc.psolve(600 * ms)
	Elapsed = run_run_run(config['tstop'])
	# time.sleep(pc.id()*2)
	
	result = {'gid_cell_id':{}, 't': list(t), 'record_soma_v':{}, 'record_spike_times':{}}
	for gid in CRC.gid_host_list:
		if pc.gid_exists(gid):
			cell_id = CRC.gid_cid[gid]
			cell = CRC.obj_cell_cid[cell_id]
			
			result['gid_cell_id'][gid] = cell_id
			# result['record_soma_v'][cell_id] = list(pc.gid2cell(gid).soma_v)

			result['record_soma_v'][cell_id] = list(cell.soma_v)
			result['record_spike_times'][cell_id] = list(cell.spike_times)
		else:
			exit(f'Error {gid} not in {CRC.gid_host_list}')

	information_print(CRC)

	with open(f'result/host_{pc.id()}.json', 'w') as outfile:
		json.dump(result, outfile)
	print(f'result/host_{pc.id()}.json 	Elapsed:{Elapsed}')

	pc.barrier()
	pc.done()
	h.quit()
