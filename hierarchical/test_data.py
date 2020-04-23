import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append('../')

from epidemics.data.combined import RegionalData
vd_data = RegionalData('VD',True)

print(len(vd_data.infected))
print(len(vd_data.time))

fig = plt.figure()
plt.plot(vd_data.infected,label='infected')
plt.plot(vd_data.recovered,label='recovered')
plt.plot(vd_data.icu,label='icu')
plt.plot(vd_data.hospitalized,label='hospitalized')
plt.plot(vd_data.deaths,label='deaths')
plt.plot(vd_data.released,label='released')
plt.plot(vd_data.ventilated,label='ventilated')
plt.legend()
plt.savefig('data.pdf')
