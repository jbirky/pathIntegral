import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
plt.style.use('classic')
rc('font', family='serif')
rc('figure', facecolor='w')
import os
import math


def phiInit(x):

	phi0 = []
	for i in range(len(x)):
		# phi0.append((2/math.pi)**(1/4) * math.exp(-2 * (x[i] - .75)**2 /2))
		phi0.append(.5*x[i]**2)

	return phi0



if __name__ == '__main__':

	# Plot expected values for each time step

	exp_files = os.listdir('expected/')

	for file in exp_files:
		vals = []
		with open('expected/' + file) as f:
			for line in f:
				vals.append(float(line))

		tstep = np.arange(0,len(vals),1)

		fname = file.split('.dat')[0]

		plt.figure(figsize=[12,6])
		plt.plot(tstep, vals)
		plt.xlabel('Timestep', fontsize=15)
		plt.ylabel(fname, fontsize=15)
		plt.savefig('plots/' + fname + '.png')
		# plt.show()
		plt.close()


	plt.figure(figsize=[12,6])
	for file in exp_files:
		vals = []
		with open('expected/' + file) as f:
			for line in f:
				vals.append(float(line))

		tstep = np.arange(0,len(vals),1)

		fname = file.split('.dat')[0]

		plt.plot(tstep, vals, label=fname, alpha=.8)
	plt.xlim(0,16)
	plt.xlabel('Timestep', fontsize=15)
	plt.ylabel('Expected Value', fontsize=15)
	plt.legend(loc='lower right')
	plt.title('Expected values vs. Time', fontsize=18)
	plt.savefig('plots/expected.png')
	# plt.show()
	plt.close()


	# Plot probability functions 

	avg_eng_array = []
	with open('expected/avg_eng.dat') as f:
		for line in f:
				avg_eng_array.append(float(line))

	avg_eng_val = np.mean(avg_eng_array)

	prob_files = os.listdir('wave_prob/')

	for file in prob_files:
		vals = []
		with open('wave_prob/' + file) as f:
			for line in f:
				vals.append(float(line))
		vals = np.array(vals)

		xvals = np.linspace(-4, 4, len(vals))

		fname = file.split('.dat')[0]
		num = fname.split('phi_sq')[1]

		plt.figure(figsize=[12,6])
		plt.plot(xvals, vals + avg_eng_val, label=r'$t_{s} = %s$'%(num), color='r')
		plt.plot(xvals, phiInit(xvals), color='k', alpha=.8, label=r'$V(x) = \frac{1}{2}m\omega^{2}$')
		plt.axhline(y=avg_eng_val, linestyle='--', color='b', alpha=.4, label=r'$<E>$')
		plt.axvline(x=.75, linestyle='--', color='k', alpha=.4)
		plt.ylim([0,2])
		plt.xlabel(r'$x$', fontsize=15)
		plt.ylabel(r'$|\Psi(x,t)|^{2}$', fontsize=15)
		plt.legend(loc='upper right')
		plt.title('Wave Probability Amplitude', fontsize=18)
		plt.savefig('plots/phi/' + fname + '.png')
		# plt.show()
		plt.close()

