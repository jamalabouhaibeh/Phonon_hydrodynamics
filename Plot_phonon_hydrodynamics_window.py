# Need to run 'Average_phonon_scattering_rates_for_multiple_temperatures.py' before using it.

import numpy as np
import pandas as pd
import matplotlib.pyplot  as plt

results = np.array(results)

plt.figure()

plt.plot(results[:,0],results[:,3], color='black', label='Normal')
plt.plot(results[:,0],results[:,4], color='red', label='Umklapp')
plt.plot(results[:,0],results[:,6], color='blue', label='Boundary at L= 500 nm')
plt.axvspan(15, 200,color='grey')


plt.xlabel("Temperature(K)")
plt.ylabel("Average Scattering Rates(1/ps)")
axes = plt.gca()
axes.set_xlim([0,300])
axes.set_ylim([1e-5,1e0])
plt.yscale('log')
plt.legend(facecolor='cyan', framealpha=1, loc=0)


plt.savefig(input("Enter the path where you would like to save the plot")+"/plot.png", dpi=300)


