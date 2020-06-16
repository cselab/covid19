#!/usr/bin/env python3

import json
import numpy as np
import argparse
import os
from datetime import datetime
from glob import glob
import re
import pandas as pd

import countrydata

# XXX path to folder with output from `request_country.py`
datafolder = "."
df = countrydata.CollectCountryData(datafolder)
countrydata.AppendOfficalLockdown(df)
countrydata.AppendInferred(df, datafolder)

with open("tint.csv", 'w') as f:
    f.write("country,inferred,official\n")
    for i, row in enumerate(df.sort_values(by='fullname').itertuples()):
        inferred = row.tint_mean.strftime('%Y-%m-%d')
        official = "" if pd.isnull(
            row.official_lockdown) else row.official_lockdown.strftime(
                '%Y-%m-%d')
        f.write("{},{},{}\n".format(row.fullname, inferred, official))
