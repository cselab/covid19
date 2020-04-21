import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from plot import Renderer

#BUILD_DIR = os.path.join(os.path.dirname(__file__), '..', 'build')
#if os.path.exists(BUILD_DIR):
#    sys.path.append(BUILD_DIR)


CANTON_TO_INDEX = {'AG': 0, 'AI': 1, 'AR': 2, 'BE': 3, 'BL': 4, 'BS': 5, 'FR': 6, 'GE': 7, 'GL': 8, 'GR': 9, 'JU': 10, 'LU': 11, 'NE': 12, 'NW': 13, 'OW': 14, 'SG': 15, 'SH': 16, 'SO': 17, 'SZ': 18, 'TG': 19, 'TI': 20, 'UR': 21, 'VD': 22, 'VS': 23, 'ZG': 24, 'ZH': 25}

NUM_CANTONS = len(CANTON_TO_INDEX)


def plot_ode_results(result,utility,n):
    """
    v = utility[number of sensors][canton][day]
    """
    days = utility.shape[2]

    def frame_callback(rend):
        t = rend.get_frame() * (days - 1) // rend.get_max_frame()
        util = utility[n,:,t]
        res  = result[t,:]
        

        v_u = {}
        v_r = {}
        texts = {}
        for i, c in enumerate(rend.get_codes()):
            i_state = CANTON_TO_INDEX[c]
            v_u[c] = util[i_state]
            v_r[c] = res [i_state]
            texts[c] = str("{:.3f}".format(v_r[c]))

        rend.set_values(v_r)
        rend.set_texts(texts)
        rend.set_bars(v_u)
        plt.suptitle("Day :" + str(t), fontsize=12)
#
    rend = Renderer(frame_callback)
    rend.save_movie(frames=days,filename=str(n) + ".mp4")


if __name__ == '__main__':
    results = np.load("output.npy")

    res = results[0,:,:]
    utility = np.load("result.npy")
    for n in range(0,utility.shape[0]):
        plot_ode_results(res,utility,n)