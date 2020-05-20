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
    "Russian Federation",
    "Ukraine",
    "France",
    "Spain",
    "Sweden",
    "Norway",
    "Germany",
    "Finland",
    "Poland",
    "Italy",
    "United Kingdom",
    "Romania",
    "Belarus",
    "Kazakhstan",
    "Greece",
    "Bulgaria",
    "Iceland",
    "Hungary",
    "Portugal",
    "Austria",
    "Czech Republic",
    "Serbia",
    "Ireland",
    "Lithuania",
    "Latvia",
    "Croatia",
    "Bosnia and Herzegovina",
    "Slovakia",
    "Estonia",
    "Denmark",
    "Switzerland",
    "Netherlands",
}

country_to_idx = {row['country']:idx for idx,row in enumerate(js)}

for country in countries:
    assert country in country_to_idx , "Error: unknown country '{}'".format(country)

for country in countries:
    idx = country_to_idx[country]
    row = js[idx]
    if country not in countries:
        continue
    data = row['confirmed']
    pop = row['population']
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
