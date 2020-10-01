#! /usr/bin/env/python3
import os
import os.path
from os import path
import csv
import json
import argparse
import numpy as np
import pandas as pd


from itertools import groupby
from operator import itemgetter


def joinres(a,b,c,d,e):
    dct0 = dict(b)
    dct1 = dict(c)
    dct2 = dict(e)
    tmp0 = dict([ (k, [v]+[dct0[k]]+list(dct1[k])) for k,v in d if k in dct1 ])
    tmp1 = dict([ (k, [v]+list(tmp0[k])) for k,v in e if k in tmp0 ])
    ret = [ (k, [v]+list(tmp1[k])) for k,v in a if k in tmp1 ]
    return ret


def findResults(base, res):
    dirs = [os.path.join(base, o) for o in os.listdir(base) 
                                if os.path.isdir(os.path.join(base,o))]

    models = [os.path.basename(d) for d in dirs]
    dirs = [os.path.join(d, res) for d in dirs]

    dirs = zip(models, dirs)
    dirs = [ (m,d) for m, d in dirs if path.exists(d)]
    res  = [ (m,os.path.join(d, 'latest')) for m, d in dirs]
    res  = [ (m,r) for m,r in res if path.exists(r)]
    return res


def getSamples(resfiles, par):
  parsamples = []
  for m, file in resfiles:
    with open(file) as f:
      r = json.load(f)
    
      idx = -1
      for i, v in enumerate(r['Variables']):
          if v['Name'] == par:
            idx = i
 
      if idx < 0:
          print("Variable not found, skip..")
          continue

      samples = np.array(r['Solver']['Posterior Sample Database'])

      if len(samples) < 1:
          print("Empty results found, skip..")
          continue

      parsamples.append( (m,samples[:,idx]) )

  return parsamples



def getStats(model, samples, pct):
    mean = np.mean(samples)
    median = np.quantile(samples,0.5)
    hi = np.quantile(samples,0.5+0.5*pct)
    lo = np.quantile(samples,0.5-0.5*pct)
    return (model, mean, median, lo, hi)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--src', type=str, help='Directory to traverse and look for result files.', required=True)
    parser.add_argument('--out', type=str, help='Output directory.', required=True)
    parser.add_argument('--res', type=str, default='_korali_samples', help='Name of sample folders.', required=False)
    parser.add_argument('--pars', type=str, default=["R0"] , nargs='+', help='List of parameter to analyse.', required=False)
    parser.add_argument('--pct', type=float, default=0.95, help='Hi/Lo Confidence interval.')
    args = parser.parse_args()

    dirs = findResults(args.src, args.res)
    print("Processing {0}..".format(args.src))
    print("{0} results found.".format(len(dirs)))

    params = args.pars
    for p in params:
        outfile="{0}/param_{1}.csv".format(args.out, p)
        print("\n\tAnalyzing param {0}.. ({1}%-CI)".format(p, args.pct))

        samples = getSamples(dirs, p)
        stats = [ getStats(m, s, args.pct) for m,s in samples ]
        df = pd.DataFrame(stats,columns=["model", "mean", "median", "low", "high"])
        print(df)
        df.to_csv(outfile, index=True, header=True)
