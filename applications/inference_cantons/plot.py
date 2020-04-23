#!/usr/bin/env python3
# Author: Petr Karnakov
# Date:   23/04/2020
# Email:  kpetr@ethz.ch

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from epidemics.tools.tools import load_model

from pathlib import Path
import argparse
import os

# required for load_model
from main import Model

dataFolder = Path("data")
f = dataFolder / 'switzerland' / 'cantons' / 'state.pickle'

assert os.path.isfile(f)

a = load_model(f)

# a.propagate()
# a.save()

a.plot_intervals()


_, samples_file = get_last_generation( args.dataFolder + "/_korali_samples/", 'gen*.json')
with open(samples_file) as json_file:
  data = json.load(json_file)
p = data["Results"]["Sample Database"];
beta  = [row[0] for row in p]
gamma = [row[1] for row in p]

fig = plt.figure(figsize=(12, 8))
ax  = fig.subplots(1)

R0 = np.asarray(beta)/np.asarray(gamma)*d['S']['Mean'][0] / d['Population Size']

ax.hist( R0 , 100, density=True, facecolor='g', alpha=0.75)
ax.set_xlabel('R0')
ax.grid()

fig.savefig(res_folder+'R0.png')

plt.show()
