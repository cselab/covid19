#!/usr/bin/env python

import json
import subprocess as sp

presets = "../../repoweb/wsgi/corona/static/countries.json"

with open(presets) as f:
    js = json.loads(f.read())

print(len(js))


def filter_few(confirmed, population):
    i = 0
    for i, c in enumerate(confirmed):
        if c > 5 and c > 2e-6 * population:
            break
    return confirmed[i:]


for row in js:
    country = row["country"]
    data = row["confirmed"]
    pop = row["population"]
    data = filter_few(data, pop)
    outdir = country.replace(' ', '')
    cmd = ["./main.py"] + \
            ["--dataFolder", outdir] + \
            ["--populationSize", str(pop)] + \
            ["--data"] + list(map(str, data))
    o = sp.run(cmd)
    print(o)
    exit()
