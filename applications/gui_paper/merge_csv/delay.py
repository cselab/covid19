#!/usr/bin/env python3

import pandas as pd

'''
For the csv file 'delay.csv'
1. Put 'merge.csv' to current directory
2. Run this script.
'''

merged = pd.read_csv('merged.csv')
out = merged.copy()

print(len(merged))

out['delay'] = len(out)*[0]

for i in range(len(out)):
  date = merged['date']
  exact = merged['exactDate']
  if not exact.empty:
      date = pd.to_datetime(date)
      exact = pd.to_datetime(exact)
      out['delay'] = (date - exact).dt.days

out.to_csv('delay.csv')
