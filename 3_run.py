from neuron import h
from neuron.units import ms, mV
from Circuit import Circuit
import time
import json

# clear; mpiexec -n 3 nrniv

num_neurons = 5

cr = Circuit(num_neurons)
pc = h.ParallelContext()
pc.set_maxstep(10 * ms)


t = h.Vector().record(h._ref_t)
h.finitialize(-65 * mV)
print(f'run №{pc.id()}')
t1 = time.time()
pc.psolve(600 * ms)
t2 = time.time()
print(f'	t({pc.id()}) = ', t2-t1)

result = {'available_gid':[], 't': list(t), 'gid':{}}
for gid in range(num_neurons):
	if pc.gid_exists(gid):
		result['available_gid'].append(gid)
		result['gid'][gid] = list(pc.gid2cell(gid).soma_v)

with open(f'result/host_{pc.id()}.json', 'w') as outfile:
	json.dump(result, outfile)
print(f'result/host_{pc.id()}.json')
pc.barrier()
pc.done()
h.quit()