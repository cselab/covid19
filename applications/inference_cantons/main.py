#!/usr/bin/env python3
# Author: Petr Karnakov
# Date:   23/04/2020
# Email:  kpetr@ethz.ch

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import argparse

import numpy as np

from ode import Ode, Sir, Seir, SeirCpp
from data import Data, get_all_canton_keys, get_data_switzerland,\
        get_data_switzerland_cantons, get_data_synthetic
from model import Model


def main():
    nSamples = 5000
    x = argparse.Namespace()
    x.dataFolder = "data/"
    x.nPropagation = 20
    x.percentages = [0.5]
    x.nThreads = 8

    #data = get_data_switzerland()
    #data = get_data_synthetic()

    keys = get_all_canton_keys()
    #keys = ['ZH', 'AG', 'TI']

    data = get_data_switzerland_cantons(keys)

    data.fit_importance = [1] * len(keys)
    key_sel = []
    #key_sel = ['ZH', 'TI']
    for k in key_sel:
        data.fit_importance[keys.index(k)] *= 100

    #ode = Sir()
    #params_to_infer = ['R0', 'gamma']

    #ode = Seir()
    ode = SeirCpp()
    #params_to_infer = []
    params_to_infer = ['R0', 'Z', 'D', 'nu', 'theta_b', 'tact']
    #params_to_infer += ['beta_corr0', 'beta_corr1', 'beta_corr2', 'beta_corr3']
    params_to_infer += ['beta_corr0', 'beta_corr1', 'beta_corr2']
    data.beta_corr_regions = {
            "beta_corr0" : ['VS', 'UR'],
            "beta_corr1" : ['AG', 'ZG', 'TG', 'JU', 'AR', 'AI', 'TI'],
            "beta_corr2" : ['BS', 'BL', 'SH', 'SO'],
            }
    data.beta_corr_regions = {k:list(map(keys.index, v)) for k,v in data.beta_corr_regions.items()}

    #ode.params_fixed['tact'] = 1e10
    #params_to_infer = ['R0', 'Z', 'D', 'nu', 'theta_b']

    #ode.params_fixed['nu'] = 0
    #params_to_infer = ['R0', 'Z', 'D', 'theta_b', 'tact']

    a = Model(data, ode, params_to_infer, **vars(x))

    if 1:
        a.sample(nSamples)
        a.propagate()
    else:
        a.evaluate([1.34, 1., 2., 2.5])

    a.save()


if __name__ == "__main__":
    main()
