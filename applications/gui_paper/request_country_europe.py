#!/usr/bin/env python

import json
import subprocess as sp
import os

presets = "./countries.json"

with open(presets) as f:
    js = json.loads(f.read())

def filter_few(confirmed, population):
    i = 0
    for i, c in enumerate(confirmed):
        if c > 5 and c > 2e-6 * population:
            break
    return confirmed[i:]


countries = {
    "Switzerland",
    "Spain",
    "Italy",
    "Austria",
    "Germany",
    "Ireland",
    "Germany",
    "Hungary",
    "Russian Federation",
    "Sweden",
    "United Kingdom",
}

for row in js:
    country = row["country"]
    if country not in countries:
        continue
    data = row["confirmed"]
    pop = row["population"]
    data = filter_few(data, pop)
    outdir = country.replace(' ', '')
    if os.path.isdir(outdir):
        print("skip existing {}".format(outdir))
        continue
    cmd = ["./main.py"] + \
            ["--dataFolder", outdir] + \
            ["--nThreads", "8"] + \
            ["--nSamples", "5000"] + \
            ["--populationSize", str(pop)] + \
            ["--data"] + list(map(str, data))
    print(cmd)
    o = sp.run(cmd)
