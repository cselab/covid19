import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append('../../')

from epidemics.data.combined import RegionalData
vd_data = RegionalData('VD',True)

fig = plt.figure()
plt.plot(vd_data.infected,label='infected')
if vd_data.recovered is not None:
    plt.plot(vd_data.recovered,label='recovered')
if vd_data.deaths is not None:
    plt.plot(vd_data.deaths,label='deaths')
if vd_data.icu is not None:
    plt.plot(vd_data.icu,label='icu')
if vd_data.hospitalized is not None:
    plt.plot(vd_data.hospitalized,label='hospitalized')
if vd_data.released is not None:
    plt.plot(vd_data.released,label='released')
if vd_data.ventilated is not None:
    plt.plot(vd_data.ventilated,label='ventilated')

plt.legend()
plt.savefig('data/data.pdf')
