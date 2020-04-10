from neuron import h
from neuron.units import ms, mV
from Circuit import Circuit
import time
import json
import numpy
import os

# clear; mpiexec -n 3 nrniv 3_run.py

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

def get_neuron_number(folder='template'):
	num_neurons = 0
	for _template_paths in os.listdir(folder):
		if _template_paths.endswith('.hoc'):
			num_neurons += 1
	return num_neurons

if __name__ == '__main__':
	config = {'tstop':600, 'iclamp_delay':50, 'iclamp_dur':100, 'iclamp_amp':0.1}
	num_neurons = get_neuron_number()

	cr = Circuit(config, num_neurons)

	pc = h.ParallelContext()
	pc.set_maxstep(10 * ms)

	t = h.Vector().record(h._ref_t)
	h.finitialize(-65 * mV)
	print(f'\nrun №{pc.id()}')

	# pc.psolve(600 * ms)
	run_run_run(config['tstop'])

	result = {'gid_cell_id':{}, 't': list(t), 'record':{}}
	for gid in range(num_neurons):
		if pc.gid_exists(gid):
			cell_id = cr.template_cell_ids[gid]
			result['gid_cell_id'][gid] = cell_id
			result['record'][cell_id] = list(pc.gid2cell(gid).soma_v)

	with open(f'result/host_{pc.id()}.json', 'w') as outfile:
		json.dump(result, outfile)
	print(f'result/host_{pc.id()}.json')
	pc.barrier()
	pc.done()
	h.quit()
