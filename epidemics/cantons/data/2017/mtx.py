#!/usr/bin/env python3

from matrixconverters.read_ptv import ReadPTVMatrix
import numpy as np
import os

p = "DWV_2017_OeV_Wegematrizen_CH_bin.mtx"
cache = os.path.splitext(p)[0] + ".npy"

if not os.path.isfile(cache):
    m = ReadPTVMatrix(filename=p)
    data = {}

    f = "zone_name"
    data[f] = m[f]

    f = "matrix"
    data[f] = m[f].astype(np.float32)

    np.save(cache, data)
    print("cache to '{:}'".format(cache))
else:
    data = np.load(cache).item()

m = data['matrix']
zz = data['zone_name']

#names = [str(z.data) for z in zz]
#coords = [int(z.coords['zone_no'].data) for z in zz]
#ii = np.argsort(names)

with open("zones", 'w') as f:
    for z in zz:
        f.write("{:} {:}\n".format(z.data, z.coords['zone_no'].data))
