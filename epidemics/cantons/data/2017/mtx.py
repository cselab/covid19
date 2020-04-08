#!/usr/bin/env python3

from matrixconverters.read_ptv import ReadPTVMatrix
import numpy as np
import cache
import geopandas as gpd


def load_matrix(path):
    cachename = path
    r = cache.load(cachename)
    if r is not None: return r

    m = ReadPTVMatrix(filename=p)
    matrix = m['matrix'].astype(np.float32)
    zonenames = [str(z.data) for z in m['zone_name']]

    r = matrix, zonenames
    return cache.save(cachename, r)

def load_zones(path):
    cachename = path
    r = cache.load(cachename)
    if r is not None: return r

    gdf = gpd.read_file(p)
    zonenames = list(map(str, gdf.N_Gem))
    cantons = list(map(str, gdf.N_KT))

    r = zonenames, cantons
    return cache.save(cachename, r)

p = "public_transport.mtx"
matrix, zonenames = load_matrix(p)

p = "zones.gpkg"
zonenames, cantons = load_zones(p)
