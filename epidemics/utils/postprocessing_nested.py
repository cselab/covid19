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


def getEvidence(resfiles):
  evidence = []
  for m, file in resfiles:
    with open(file) as f:
      r = json.load(f)
    
      e = r['Solver']['LogEvidence']
      evidence.append( (m, e) )

  return evidence

def getVariance(resfiles):
  variance = []
  for m, file in resfiles:
    with open(file) as f:
      r = json.load(f)
    
      e = r['Solver']['LogEvidence Var']
      variance.append( (m, e) )

  return variance


def computeProbabilities(evidence):
    maxlogevidence = -999999999999
    for m, e in evidence:
        if e > maxlogevidence:
            maxlogevidence = e

    evidence = [(m, np.exp(e - maxlogevidence)) for m, e in evidence]
 
    sumevidence = 0.0
    for m, e in evidence:
        sumevidence += e

    ps = [(m, e/sumevidence) for m, e in evidence]
    return ps


def getBestLLk(resfiles):
  best = []
  for m, file in resfiles:
    with open(file) as f:
      r = json.load(f)
    
      e = r['Solver']['Max Evaluation']
      best.append( (m, e) )

  return best


def getStats(samples, pct):
    mean = np.mean(samples)
    median = np.quantile(samples,0.5)
    hi = np.quantile(samples,0.5+0.5*pct)
    lo = np.quantile(samples,0.5-0.5*pct)
    return (mean, median, hi, lo)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--src', type=str, help='Directory to traverse and look for result files.', required=True)
    parser.add_argument('--res', type=str, default='_korali_samples', help='Name of sample folders.', required=False)
    parser.add_argument('--par', type=str, help='Name of parameter to analyse.', required=True)
    parser.add_argument('--out', type=str, default='nested_post.csv', help='Output file.')
    args = parser.parse_args()

    dirs = findResults(args.src, args.res)
    print("Processing {0}..".format(args.src))
    print("{0} results found.".format(len(dirs)))

    samples = getSamples(dirs, args.par)
    stats = [ (m,getStats(s, 0.9)) for m,s in samples]

    evidence = getEvidence(dirs)
    probs = computeProbabilities(evidence)
    variance = getVariance(dirs)
    best = getBestLLk(dirs)

    out = joinres(evidence, best, stats, probs, variance)
    
    df = pd.DataFrame(out,columns=["model", "values"])
    df3 = pd.DataFrame(df["values"].tolist(), columns=['evidence', 'variance', 'prob', 'best', 'mean', 'median', 'high', 'low'], index=df.model)
    print(df3)
    df3.to_csv(args.out, index=True, header=True)
