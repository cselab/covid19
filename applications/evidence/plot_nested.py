#!/usr/bin/env python3

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

import argparse
import pickle
from epidemics.tools.nested import plotNetsedResult


if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--resultFile', '-rf', required=True, type=str, help='Datafile for synthetic data.')


    args = parser.parse_args()
    filename = args.resultFile
    result   = pickle.load( open(filename, "rb") )
    
    plotNetsedResult(result)
