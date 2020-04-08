#!/usr/bin/env python3

from mpl_toolkits.basemap import Basemap
import numpy as np

# Zurich
lat0 = 47.3769; lon0 = 8.5417
m = Basemap(projection='aeqd',
            lon_0 = lon0,
            lat_0 = lat0,
            width = 1,
            height = 1)

shpfile = "cantons"
shp = m.readshapefile(shpfile, 'shapes')

d = {}
for info, shape in zip(m.shapes_info, m.shapes):
    x, y = zip(*shape)
    xy = np.array((x, y)).astype(np.float32)
    name = info['NAME']
    if name not in d:
        d[name] = []
    d[name].append(xy)

np.save("canton_shapes.npy", d)

