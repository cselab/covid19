#!/usr/bin/env python3

"""Plot reference data."""

import os
import sys
import numpy as np

from data import fetch_canton_data
from plot import Renderer

CANTON_TO_INDEX, IR = fetch_canton_data()
IR_MAX = np.nanmax(IR)

def plot_data():
    def frame_callback(rend):
        values = rend.get_values()
        texts = rend.get_texts()
        frame = rend.get_frame()
        for i, c in enumerate(rend.get_codes()):
            Ir = IR[frame][CANTON_TO_INDEX[c]]
            if np.isfinite(Ir):
                values[c] = Ir / IR_MAX * 2
                texts[c] = str(int(Ir))
        rend.set_values(values)
        rend.set_texts(texts)

    rend = Renderer(frame_callback)
    rend.save_movie(frames=len(IR), filename="data.mp4", fps=5)
    rend.save_image(filename="data.png")


def main():
    plot_data()

if __name__ == '__main__':
    main()
