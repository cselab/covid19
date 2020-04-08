#!/usr/bin/env python3
# Author: Ivica Kicci
# Date:   2020-04-06
# Email:  kicici@ethz.ch

print("DEPRECATED!")

import io
import math
import numpy as np
import os
import pandas as pd
import requests
import sys
import random
from scipy.integrate import solve_ivp

BUILD_DIR = os.path.join(os.path.dirname(__file__), '..', 'build')
if os.path.exists(BUILD_DIR):
    sys.path.append(BUILD_DIR)

from  .data import fetch_canton_data, CANTON_POPULATION, get_symmetric_Mij
import libsolver

from ..std_models.std_models import standard_deviation_models, standardDeviationModelConst
from ..epidemics import epidemicsBase
from ..tools.tools import save_file
from .misc import Values, flatten, filter_out_nans_wrt, flatten_and_remove_nans


class MultiRegionModel( epidemicsBase ):
    def __init__( self, fileName=[], defaultProperties={}, **kwargs ):
        # Ignore the following arguments:
        kwargs.pop('country')
        kwargs.pop('rawData')
        kwargs.pop('populationSize')
        kwargs.pop('nPropagation')
        kwargs.pop('stdModel')

        self.modelName        = 'multi_seiin'
        self.modelDescription = "Multi-region SEII model."

        defaultProperties = { **defaultProperties,
            'futureDays': 10,
            'dataFolder': './data/',
            'nValidation': 0
        }

        super().__init__( fileName=fileName, defaultProperties=defaultProperties, **kwargs )

        if fileName == []:
            self.download_raw_data()
            self.propagationData={}
            self.process_data()

    def download_raw_data( self ):
        cantons, infected = fetch_canton_data()
        num_days = len(infected)
        self.numCantons = len(cantons)
        self.data['Raw']['Cantons'] = cantons
        self.data['Raw']['Time'] = np.asarray(range(num_days))
        self.data['Raw']['Infected'] = infected

        Mij = get_symmetric_Mij(cantons)
        assert len(Mij) == len(cantons), (len(Mij), len(cantons))
        assert len(Mij[0]) == len(cantons)
        self.data['Raw']['Flat Mij'] = flatten(Mij)

    def set_variables_and_distributions( self ):
        p = [ 'beta', 'mu', 'alpha', 'Z', 'D', 'theta', '[Sigma]' ]

        for k, x in enumerate(p):
            self.e['Variables'][k]['Name'] = x
            self.e['Variables'][k]['Prior Distribution'] = 'Prior for ' + x

        self.nParameters = len(p)

        # Note: the order here must match the order in src/model.h:Parameters!
        k = 0
        self.e['Distributions'][k]['Name'] = 'Prior for beta'
        self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
        self.e['Distributions'][k]['Minimum'] = 0.8
        self.e['Distributions'][k]['Maximum'] = 1.5
        k += 1

        self.e['Distributions'][k]['Name'] = 'Prior for mu'
        self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
        self.e['Distributions'][k]['Minimum'] = 0.2
        self.e['Distributions'][k]['Maximum'] = 1.0
        k += 1

        self.e['Distributions'][k]['Name'] = 'Prior for alpha'
        self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
        self.e['Distributions'][k]['Minimum'] = 0.02
        self.e['Distributions'][k]['Maximum'] = 1.0
        k += 1

        self.e['Distributions'][k]['Name'] = 'Prior for Z'
        self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
        self.e['Distributions'][k]['Minimum'] = 2.0  # Days.
        self.e['Distributions'][k]['Maximum'] = 5.0  # Days.
        k += 1

        self.e['Distributions'][k]['Name'] = 'Prior for D'
        self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
        self.e['Distributions'][k]['Minimum'] = 2.0  # Days.
        # self.e['Distributions'][k]['Maximum'] = 5.0  # Days.
        self.e['Distributions'][k]['Maximum'] = 15.0  # Days.
        k += 1

        self.e['Distributions'][k]['Name'] = 'Prior for theta'
        self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
        self.e['Distributions'][k]['Minimum'] = 1.0
        self.e['Distributions'][k]['Maximum'] = 1.75
        k += 1

        self.e['Distributions'][k]['Name'] = 'Prior for [Sigma]'
        self.e['Distributions'][k]['Type'] = 'Univariate/Uniform'
        self.e['Distributions'][k]['Minimum'] = 0.0
        self.e['Distributions'][k]['Maximum'] = +600.0

    def save_data_path( self ):
        return ( self.dataFolder, self.modelName )

    def process_data( self ):
        y = self.data['Raw']['Infected']
        t = self.data['Raw']['Time']
        cantons = self.data['Raw']['Cantons']
        numCantons = len(cantons)

        N0 = [CANTON_POPULATION[canton] for canton in cantons]
        S0 = N0
        E0 = [0] * numCantons
        IR0 = [0] * numCantons
        IU0 = [0] * numCantons

        # Ticino.
        IR0[cantons['TI']] = 1
        IU0[cantons['TI']] = 10

        y0 = S0 + E0 + IR0 + IU0 + N0

        # Only non-nan data is used for evaluation.
        numDays = len(t) - self.nValidation
        self.data['Model']['x-data'] = t[1:numDays]
        self.data['Model']['y-data'] = flatten_and_remove_nans(y[1:numDays])
        self.data['Model']['Full y-data'] = flatten(y[1:numDays])
        if self.nValidation > 0:
            self.data['Validation']['x-data'] = t[-self.nValidation:]
            self.data['Validation']['y-data'] = flatten_and_remove_nans(y[-self.nValidation:])

        self.data['Model']['Initial Condition'] = y0

        T = t[-1] + self.futureDays
        self.data['Propagation']['x-data'] = list(range(T))

        self.solver = libsolver.Solver(self.data['Raw']['Flat Mij'])
        save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )

    def computational_model( self, s ):
        p = s['Parameters']
        t  = self.data['Model']['x-data']
        y0 = self.data['Model']['Initial Condition']

        # beta, mu, alpha, Z, D, theta, [sigma]
        params = libsolver.Parameters(*p[:-1])
        result_all = self.solver.solve(params, y0, int(t[-1]))
        result_Ir = [states.Ir() for state in result_all]

        y, d = filter_out_nans_wrt(flatten(result_Ir), self.data['Model']['Full y-data'])
        s['Reference Evaluations'] = y
        s['Standard Deviation Model'] = [p[-1] for known_data_cell in y]


    def computational_model_propagate( self, s ):
        p = s['Parameters']
        t  = self.data['Propagation']['x-data']
        y0 = self.data['Model']['Initial Condition']

        params = libsolver.Parameters(*p[:-1])
        result_all = self.solver.solve(params, y0, int(t[-1]))
        result_S  = [state.S()  for state in result_all]
        result_Ir = [state.Ir() for state in result_all]

        js = {}
        js['Variables'] = [{}, {}]
        js['Variables'][0]['Name'] = 'S'
        js['Variables'][0]['Values'] = flatten(result_S)
        js['Variables'][1]['Name'] = 'I'
        js['Variables'][1]['Values'] = flatten(result_Ir)

        js['Number of Variables'] = len(js['Variables'])
        js['Length of Variables'] = len(js['Variables'][0]['Values'])

        y = flatten(result_Ir[:-self.futureDays])
        js['Standard Deviation Model'] = [p[-1] for known_data_cell in y]
        s['Saved Results'] = js

    def set_variables_for_interval( self ):
        self.intervalVariables = {}

    def plot_intervals( self ):
        # TODO: Use the plotting script.
        return
