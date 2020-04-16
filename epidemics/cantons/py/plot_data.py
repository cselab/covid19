#!/usr/bin/env python3

"""Plot reference data."""

import numpy as np

import json
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from epidemics.cantons.py.model import get_canton_model_data, get_canton_reference_data
from epidemics.cantons.py.plot import Renderer
from epidemics.data.swiss_cantons import json_to_numpy_matrix
from epidemics.tools.tools import flatten

IR = get_canton_reference_data().cases_per_country
IR_MAX = np.nanmax(flatten(IR.values()))

# TODO: Move to epidemics/data/**
DATA_DIR = os.path.join(os.path.normpath(os.path.dirname(__file__)), '..', 'data')

def plot_data():
    def frame_callback(rend):
        values = rend.get_values()
        texts = rend.get_texts()
        frame = rend.get_frame()
        for i, c in enumerate(rend.get_codes()):
            Ir = IR[c][frame]
            if np.isfinite(Ir):
                values[c] = Ir / IR_MAX * 2
                texts[c] = str(int(Ir))
        rend.set_values(values)
        rend.set_texts(texts)

    # First render the default model data (both Mij and Cij).
    model_data = get_canton_model_data()
    rend = Renderer(frame_callback, data=model_data)
    rend.save_image(filename="data.png")
    rend.save_movie(frames=len(next(iter(IR.values()))), filename="data.mp4", fps=5)

    # Then the data from the other source (only Mij).
    with open(os.path.join(DATA_DIR, '2017', 'matrix.json')) as f:
        Mij_json = json.load(f)
    del Mij_json['Enk.']
    del Mij_json['LIE']
    model_data.Mij = json_to_numpy_matrix(Mij_json, model_data.region_keys)
    model_data.Cij *= 0.0
    rend = Renderer(frame_callback, data=model_data)
    rend.save_image(filename="data_2017.png")


def main():
    plot_data()

if __name__ == '__main__':
    main()
