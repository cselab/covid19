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


def joinres(a,b,c):
    dct0 = dict(b)
    dct1 = dict(c)
    ret = [ (k, [v]+[dct0[k]]+list(dct1[k])) for k,v in a if k in dct1 ]
    return ret


def findResults(base):
    dirs = [os.path.join(base, o) for o in os.listdir(base) 
                                if os.path.isdir(os.path.join(base,o))]

    models = [os.path.basename(d) for d in dirs]

    dirs = zip(models, dirs)
    dirs = [ (m,d) for m, d in dirs if path.exists(d)]
    res  = [ (m,os.path.join(d, 'cmaes.json')) for m, d in dirs]
    res  = [ (m,r) for m,r in res if path.exists(r)]
    return res


def getBestLLk(resfiles):
  best = []
  for m, file in resfiles:
    with open(file) as f:
      r = json.load(f)
    
      llk = r['Value']
      p   = r['Parameter']
      best.append( (m, llk, p) )

  return best



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--src', type=str, help='Directory to traverse and look for cmaes.json result files.', required=True)
    parser.add_argument('--out', type=str, default='cmaes_post.csv', help='Output file.')
    args = parser.parse_args()

    dirs = findResults(args.src)
    print("{0} results found.".format(len(dirs)))

    best = getBestLLk(dirs)
    df = pd.DataFrame(best,columns=["model", "best", "params"])
    df = df.set_index("model")
    
    print(df)
    df.to_csv(args.out, index=True, header=True)
