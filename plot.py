import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
plt.style.use('classic')
rc('font', family='serif')
rc('figure', facecolor='w')
import os


if __name__ == '__main__':

	out_files = os.listdir('output')

	for file in out_files:
		vals = []
		with open('output/' + file) as f:
			for line in f:
				vals.append(line)

		tstep = np.arange(0,len(vals),1)

		plt.figure(figsize=[8,6])
		plt.plot(tstep, vals)
		plt.xlabel('Timestep', fontsize=15)
		plt.ylabel(file.split('.dat')[0], fontsize=15)
		plt.show()
		plt.close()