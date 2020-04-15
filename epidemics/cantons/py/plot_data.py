#!/usr/bin/env python3

"""Plot reference data."""

import numpy as np

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

from epidemics.cantons.py.model import get_canton_reference_data
from epidemics.cantons.py.plot import Renderer
from epidemics.tools.tools import flatten

IR = get_canton_reference_data().cases_per_country
IR_MAX = np.nanmax(flatten(IR.values()))

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

    rend = Renderer(frame_callback)
    rend.save_image(filename="data.png")
    rend.save_movie(frames=len(next(iter(IR.values()))), filename="data.mp4", fps=5)

    rend = Renderer(frame_callback, matrix_json='2017/matrix.json')
    rend.save_image(filename="data_2017.png")


def main():
    plot_data()

if __name__ == '__main__':
    main()
