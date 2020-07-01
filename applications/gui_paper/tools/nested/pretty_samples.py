#!/usr/bin/env python3

from glob import glob
import os
import argparse


pretty_dict = {
    'I0': r'$k_I$',
    'R0': r'$R^\\star_0$',
    '[Dispersion]': r'[Dispersion]',
    'gamma': r'$\\gamma$',
    'tint': r'$t_{\\mathrm{int}}$',
    'dint': r'$\\delta_{\\mathrm{int}}$',
    'kbeta': r'$k_{\\mathrm{int}}$',
    'tint2_minus_tint': r'$t^{(2)}_{\\mathrm{int}}-t_{\\mathrm{int}}$',
    'dint2': r'$\\delta^{(2)}_{\\mathrm{int}}$',
    'kbeta2_div_kbeta': r'$k^{(2)}_{\\mathrm{int}}/k_{\\mathrm{int}}$',
}

parser = argparse.ArgumentParser(
    description=
    "Creates a copy of '_korali_samples' replacing the  names of variables using 'sir.prertty_dict'"
)
aa = parser.add_argument
aa('samples_dir', default='data', help="Path to '_korali_samples/'")
aa('new_samples_dir',
   default='data',
   help="Path to output (e.g. '_korali_samples_pretty'")
args = parser.parse_args()

os.makedirs(args.new_samples_dir, exist_ok=True)
assert os.path.isdir(args.samples_dir)

for fpath in glob(os.path.join(args.samples_dir, "gen*.json")):
    with open(fpath, 'r') as f:
        text = f.read()
    for varname, pretty in pretty_dict.items():
        pattern = '"Name": "{}"'
        text = text.replace(pattern.format(varname), pattern.format(pretty))
    fbase = os.path.basename(fpath)
    foutpath = os.path.join(args.new_samples_dir, fbase)
    print(foutpath)
    with open(foutpath, 'w') as fout:
        fout.write(text)
