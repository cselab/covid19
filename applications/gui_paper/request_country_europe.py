#!/usr/bin/env python

import json
import subprocess as sp

presets = "./countries.json"

with open(presets) as f:
    js = json.loads(f.read())

print(len(js))


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
}

for row in js:
    country = row["country"]
    if country not in countries:
        continue
    data = row["confirmed"]
    pop = row["population"]
    data = filter_few(data, pop)
    outdir = country.replace(' ', '')
    cmd = ["./main.py"] + \
            ["--dataFolder", outdir] + \
            ["--nThreads", "8"] + \
            ["--populationSize", str(pop)] + \
            ["--data"] + list(map(str, data))
    o = sp.run(cmd)
    print(cmd)
    #print(o)
    #exit()
