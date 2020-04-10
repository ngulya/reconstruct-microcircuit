import matplotlib.pyplot as plt
import json
import os

num_cells = 5

t = []
l_soma_v = [0 for i in range(num_cells)]

for file in os.listdir('result'):
	if file.endswith('.json'):
		with open(f'result/{file}', 'r') as outfile:
			result = json.load(outfile)
		if not len(t):
			t = result['t']
		for gid, soma_v in result['gid'].items():
			gid = int(gid)
			print(l_soma_v[gid])
			l_soma_v[gid] = soma_v

for i in range(num_cells):
	label = str(i)
	if i == 3:
		label = 'Main: ' + label
	plt.plot(t, l_soma_v[i], label = label)

plt.legend()
plt.show()
