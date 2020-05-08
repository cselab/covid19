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
    nSamples = 1000
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
    #for k in keys:
    #    i = keys.index(k)
    #    data.fit_importance[i] *= data.total_infected[i][-1] / 1000.

    #ode = Sir()
    #params_to_infer = ['R0', 'gamma']

    #ode = Seir()
    ode = SeirCpp()
    params_to_infer = []
    #params_to_infer += ['R0']
    #params_to_infer = ['R0', 'Z', 'D', 'nu', 'theta_b', 'tact']

    if 1:
        #params_to_infer += ['beta_corr' + str(i) for i in range(5)]
        data.beta_corr_regions = {
            "beta_corr0" : ['TI', 'JU'],
            "beta_corr1" : ['AG', 'SO', 'TG', 'AR'],
            "beta_corr2" : ['VS', 'UR'],
            "beta_corr3" : ['BL', 'BS', 'SH'],
            "beta_corr4" : ['ZG', 'TG'],
        }
        params_to_infer += data.beta_corr_regions.keys()

    if 0:
        data.beta_corr_regions = dict()
        keyssort = list(
            sorted(keys, key=lambda k: -data.total_infected[keys.index(k)][-1]))
        for i, k in enumerate(keyssort[:13]):
            var = "beta_corr{:}".format(i)
            data.beta_corr_regions[var] = [k]
            params_to_infer += [var]

    data.beta_corr_regions = {
        k: list(map(keys.index, v))
        for k, v in data.beta_corr_regions.items()
    }

    a = Model(data, ode, params_to_infer, **vars(x))

    if 1:
        a.sample(nSamples)
        a.propagate()
    else:
        a.evaluate([1.34, 1., 2., 2.5])

    a.save()


if __name__ == "__main__":
    main()
