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


def joinres(a,b):
    dct = dict(b)
    ret = [ (k, [v]+list(dct[k])) for k,v in a if k in dct]
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


def getEvidence(resfiles):
  evidence = []
  for m, file in resfiles:
    with open(file) as f:
      r = json.load(f)
    
      e = r['Solver']['LogEvidence']
      evidence.append( (m, e) )

  return evidence



def getStats(samples, pct):
    mean = np.mean(samples)
    median = np.quantile(samples,0.5)
    hi = np.quantile(samples,0.5+0.5*pct)
    lo = np.quantile(samples,0.5-0.5*pct)
    return (mean, median, hi, lo)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--src', type=str, help='Directory to traverse and look for result files.', required=True)
    parser.add_argument('--res', type=str, help='Name of sample folders.', required=True)
    parser.add_argument('--par', type=str, help='Name of parameter to analyse.', required=True)
    parser.add_argument('--out', type=str, default='output.csv', help='Output file.')
    args = parser.parse_args()

    dirs = findResults(args.src, args.res)
    print("{0} results found.".format(len(dirs)))

    samples = getSamples(dirs, args.par)
    stats = [ (m,getStats(s, 0.9)) for m,s in samples]
#    print(stats)

    evidence = getEvidence(dirs)
#    print(evidence)

    out = joinres(evidence, stats)
    print(out)
    
    df = pd.DataFrame(out,columns=["model", "values"])
    df3 = pd.DataFrame(df["values"].tolist(), columns=['evidence', 'mean', 'median', 'high', 'low'], index=df.model)
    print(df3)
    df3.to_csv(args.out, index=True, header=True)
