#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   16/3/2020
# Email:  garampat@ethz.ch

import glob
import importlib
import json
import os
import pickle
import shutil
import sys
from scipy.stats import truncnorm

LOGPREFIX = '[Epidemics] '

def printlog(msg, prefix=LOGPREFIX, end='\n', flush=False):
    out = sys.stdout
    out.write(f"{prefix}{msg}{end}")
    if flush:
        out.flush()

def abort(msg, code=1, prefix=LOGPREFIX):
    out = sys.stdout
    out.write(f"{LOGPREFIX}{msg}")
    out.flush()
    exit(code)

def flatten(matrix):
    """
    >>> flatten([[10, 20, 30], [40, 50]])
    [10, 20, 30, 40, 50]
    """
    return [
        value
        for row in matrix
        for value in row
    ]


def prepare_folder( dir, clean=False ):
  dir = os.path.relpath( dir )
  if( os.path.commonpath( ['data',dir])  != 'data' ):
     sys.exit(f'\n[prepare_folder] Error: The data must be saved in the ./data folder. Chosen folder: {dir}\n')
  # XXX dangerous!
  if(clean==True):
    shutil.rmtree(dir, ignore_errors=True)
  if os.path.isdir(dir) == False:
    os.makedirs(dir)




def get_last_generation( folder, pattern ):
  files = glob.glob(folder + pattern)
  return len(files), sorted(files)[-1]




def save_file( data, file, str, fileType='pickle' ):
  prepare_folder( os.path.dirname(file), clean=False )

  if str:
    print(f'[Epidemics] Save {str} in {fileType} file: {file}')

  if(fileType=='json'):
    with open(file, 'w') as f:
      json.dump( data, f, indent=2, sort_keys=True)

  if(fileType=='pickle'):
    with open(file, 'wb') as f:
      pickle.dump( data, f )




def load_file( file, str, fileType='pickle' ):

  if str:
    print(f'[Epidemics] Load {str} from {fileType} file: {file}')

  if(fileType=='json'):
    with open(file, 'r') as f:
      data = json.load( f )

  if(fileType=='pickle'):
    with open(file, 'br') as f:
      data = pickle.load( f )

  return data


def load_model(path):
    """Unpickle a model instance from the given path."""
    with open(path, 'rb') as f:
        model = pickle.load(f)
    return model



def make_path( path, *paths):
  p = os.path.join(path,*paths);
  p = os.path.normpath(p);
  return p



def import_from( module, name ):
  module = importlib.import_module( module )
  return getattr( module, name )



def get_truncated_normal(mean, sd, low, upp):
    return truncnorm( (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
