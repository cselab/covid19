#!/usr/bin/env python3

import geopandas as gpd

p = "Verkehrszonen_Schweiz_NPVM_2017.gpkg"

gdf = gpd.read_file(p)

codes = gdf.N_KT
names = gdf.N_Gem

print(gdf.keys())

with open("codes", 'w') as f:
    for name,code in zip(names,codes):
        f.write("{:} {:}\n".format(name, code))
