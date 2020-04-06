#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   31/3/2020
# Email:  garampat@ethz.ch

import io
import math
import numpy as np
import os
import pandas as pd
import requests
import sys
import random
from scipy.integrate import solve_ivp

BUILD_DIR = os.path.join(os.path.dirname(__file__), 'build')
if os.path.exists(BUILD_DIR):
    sys.path.append(BUILD_DIR)

from  .data import fetch_canton_data, CANTON_POPULATION, get_symmetric_Mij
import libseiir

from ..std_models.std_models import standard_deviation_models, standardDeviationModelConst
from ..epidemics import epidemicsBase
from ..tools.tools import save_file
from .misc import Values, flatten, filter_out_nans_wrt, flatten_and_remove_nans, extract_values_from_state


class MultiSEIIRModel( epidemicsBase ):
    def __init__( self, fileName=[], defaultProperties={}, **kwargs ):
        # Ignore the following arguments:
        kwargs.pop('country')
        kwargs.pop('rawData')
        kwargs.pop('populationSize')

        self.modelName        = 'multi_seiir'
        self.modelDescription = "Multi-region SEIIR model."

        defaultProperties = { **defaultProperties,
            'stdModel': 0,
            'futureDays': 10,
            'nPropagation': 100,
            'logPlot': False,
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
        self.data['Model']['Standard Deviation Model'] = self.stdModel

        T = t[-1] + self.futureDays
        self.data['Propagation']['x-data'] = np.linspace(0, T, self.nPropagation).tolist()

        self.multiseiir = libseiir.MultiSEIIR(self.data['Raw']['Flat Mij'])
        save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )

    def computational_model( self, s ):
        p = s['Parameters']
        t  = self.data['Model']['x-data']
        y0 = self.data['Model']['Initial Condition']

        # beta, mu, alpha, Z, D, theta, [sigma]
        params = libseiir.Parameters(*p[:-1])
        result_all = self.multiseiir.solve(params, y0, int(t[-1]))
        result_Ir = [extract_values_from_state(state, self.numCantons, Values.Ir) for state in result_all]

        y, d = filter_out_nans_wrt(flatten(result_Ir), self.data['Model']['Full y-data'])
        s['Reference Evaluations'] = y
        s['Standard Deviation Model'] = [p[-1] for known_data_cell in y]


    def computational_model_propagate( self, s ):
        p = s['Parameters']
        t  = self.data['Propagation']['x-data']
        y0 = self.data['Model']['Initial Condition']

        params = libseiir.Parameters(*p[:-1])
        result_all = self.multiseiir.solve(params, y0, int(t[-1]))
        result_S  = [extract_values_from_state(state, self.numCantons, Values.S)  for state in result_all]
        result_Ir = [extract_values_from_state(state, self.numCantons, Values.Ir) for state in result_all]

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

        self.intervalVariables['Total Infected'] = {}
        self.intervalVariables['Total Infected']['Formula'] = lambda v: self.populationSize - v['S']

        self.intervalVariables['Infected Rate'] = {}
        self.intervalVariables['Infected Rate']['Formula'] = lambda v: self.parameters[0]['Values'] * v['S'] * v['I'] / self.populationSize




    def plot_intervals( self ):
        fig = plt.figure(figsize=(12, 8))
        fig.suptitle(self.modelDescription)
        ax  = fig.subplots(len(self.credibleIntervals))
        ax[0].plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label="Total Infected (data)", color='black')

        if self.nValidation > 0:
            ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label="Total Infected (validation data)", color='black')

        z = np.asarray(self.data['Model']['y-data'])[1:] - np.asarray(self.data['Model']['y-data'])[0:-1]
        ax[1].plot( self.data['Model']['x-data'][1:], z, 'o', lw=2, label="Infected Rate (data)", color='black')

        for k,y in enumerate( self.credibleIntervals.keys() ):
            ax[k].plot( self.data['Propagation']['x-data'], self.credibleIntervals[y]['Mean'],   '-', lw=2, label="Mean", color='blue' )
            ax[k].plot( self.data['Propagation']['x-data'], self.credibleIntervals[y]['Median'], '-', lw=2, label="Median", color='black')

            self.credibleIntervals[y]['Intervals'].sort(key = lambda x: x['Percentage'], reverse = True)

            for x in self.credibleIntervals[y]['Intervals']:
                p1 = [ max(k,0) for k in x['Low Interval'] ]
                p2 = x['High Interval']
                p  = 100.*x['Percentage']
                ax[k].fill_between( self.data['Propagation']['x-data'], p1 , p2,  alpha=0.5, label=f" {p:.1f}% credible interval" )

            ax[k].legend(loc='upper left')
            ax[k].set_ylabel( y )
            ax[k].set_xticks( range( np.ceil( max( self.data['Propagation']['x-data'] )+1 ).astype(int) ) )
            ax[k].grid()
            if self.logPlot:
                ax[k].set_yscale('log')

        ax[-1].set_xlabel("time in days")

        file = os.path.join(self.saveInfo['figures'],'prediction.png');
        prepare_folder( os.path.dirname(file) )
        fig.savefig(file)

        plt.show()

        plt.close(fig)


        fig = plt.figure(figsize=(12, 8))
        ax  = fig.subplots(1)

        R0 = self.parameters[0]['Values'] / self.parameters[1]['Values']

        ax.hist( R0 , 100, density=True, facecolor='g', alpha=0.75)
        ax.set_xlabel("R0")
        ax.grid()

        file = os.path.join(self.saveInfo['figures'], 'R0.png')
        prepare_folder( os.path.dirname(file), clean=False )
        fig.savefig(file)
