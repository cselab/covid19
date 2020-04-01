#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   16/3/2020
# Email:  garampat@ethz.ch

import os
import sys
import shutil
import glob
import importlib
import pickle # XXX use cpickle
import json



def prepare_folder( dir, clean=True ):
  dir = os.path.relpath( dir )
  if( os.path.commonpath( ['data',dir])  != 'data' ):
     sys.exit('Error: The data must be saved in the ./data folder')
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

  if not str:
    print(f'\nSave {str} in  {fileType} file: {file}\n')

  if(fileType=='json'):
    with open(file, 'w') as f:
      json.dump( data, f, indent=2, sort_keys=True)

  if(fileType=='pickle'):
    with open(file, 'wb') as f:
      pickle.dump( data, f )




def load_file( file, str, fileType='pickle' ):

  if not str:
    print(f'Load {str} from {fileType} file: {file}')

  if(fileType=='json'):
    with open(file, 'r') as f:
      data = json.load( f )

  if(fileType=='pickle'):
    with open(file, 'br') as f:
      data = pickle.load( f )

  return data




def make_path( path, *paths):
  p = os.path.join(path,*paths);
  p = os.path.normpath(p);
  return p



def import_from( module, name ):
  module = importlib.import_module( module )
  return getattr( module, name )
