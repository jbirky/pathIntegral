import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
plt.style.use('classic')
rc('font', family='serif')
rc('figure', facecolor='w')
import os


if __name__ == '__main__':

	# Plot expected values for each time step

	exp_files = os.listdir('expected/')

	for file in exp_files:
		vals = []
		with open('expected/' + file) as f:
			for line in f:
				vals.append(line)

		tstep = np.arange(0,len(vals),1)

		fname = file.split('.dat')[0]

		plt.figure(figsize=[12,6])
		plt.plot(tstep, vals)
		plt.xlabel('Timestep', fontsize=15)
		plt.ylabel(fname, fontsize=15)
		plt.show()
		plt.savefig('plots/' + fname + '.png')
		plt.close()


	# Plot probability functions 

	prob_files = os.listdir('wave_prob/')

	for file in prob_files:
		vals = []
		with open('wave_prob/' + file) as f:
			for line in f:
				vals.append(line)

		xvals = np.linspace(-4, 4, len(vals))

		fname = file.split('.dat')[0]

		plt.figure(figsize=[12,6])
		plt.plot(xvals, vals)
		plt.xlabel('Timestep', fontsize=15)
		plt.ylabel(fname, fontsize=15)
		plt.show()
		plt.savefig('plots/' + fname + '.png')
		plt.close()

