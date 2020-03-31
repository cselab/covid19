#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch


import argparse

from epidemics.tools.tools import import_from

parser = argparse.ArgumentParser()
parser.add_argument('--compModel', '-cm', default='sir', help='The computational mode.')
parser.add_argument('--dataFolder', '-df', default='data/', help='Save all results in the folder \'data\\dataFolder\' ')
parser.add_argument('--country', '-c', default='switzerland', help='Country from which to retrieve data./')
parser.add_argument('--rawData', '-d', default=[], nargs='+', type=float, help='Infected population.')
parser.add_argument('--populationSize', '-ps', type=int, default=80000, help='Total population.')
parser.add_argument('--stdModel', '-sm', type=int, default=0, help='Standard deviation model. 0-Constant, 1-Sqrt, 2-Linear')
parser.add_argument('--nSamples', '-ns', type=int, default=2000, help='Number of samples for TMCMC.')
parser.add_argument('--nThreads', '-nt', type=int, default=1, help='Number of threads.')
parser.add_argument('--nPropagation', '-np', type=int, default=100, help='Number of points to evaluate the solution in the propagation phase.')
parser.add_argument('--futureDays', '-fd', type=int, default=2, help='Propagate that many days in future, after the time of observation of the last data.')
parser.add_argument('--nValidation', '-nv', type=int, default=0, help='Use that many data from the end of the data list to validate the prediction.')
parser.add_argument('--percentages', '-p', nargs='+', type=float, default=[0.5, 0.95], help='Percentages for confidence intervals.')
parser.add_argument('--silent', action='store_true', help='No output on screen.')
args = parser.parse_args()


model = import_from( 'epidemics.' + args.compModel, 'epModel')


a = model( **vars(args) )

a.sample( args.nSamples )

a.propagate()

a.compute_intervals()

a.plot_intervals()

a.save()
