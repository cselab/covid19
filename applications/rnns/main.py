#!/usr/bin/env python

import json
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
from model import Model
from data import get_data

if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--country', '-c', default='Test')
    parser.add_argument('--n_features', '-nf', default=1,type=int)
    parser.add_argument('--seq_length', '-sl', default=20,type=int)
    parser.add_argument('--hidden_dim', '-hd', default=10,type=int)
    parser.add_argument('--n_layers', '-nl', default=1,type=int)
    parser.add_argument('--n_epochs', '-ne', default=20000,type=int)
    parser.add_argument('--lr', '-lr', default=0.001,type=int)
    parser.add_argument('--name', '-n', default='',type=str)

    params = vars(parser.parse_args())

    data_source = 'synthetic'

    # Synthetic
    r0 = 1.5
    gamma = 0.1
    tact = 40
    dtact = 10
    kbeta = 0.4
    T = 100
    p = [r0, gamma, tact, dtact, kbeta]
    train_x,train_y = generate_data('synthetic',p,T)

    # # JSON
    # train_x, train_y = get_data(params['country'])
    # # Pickle
    # path = './data/SIR_single.pickle'
    # train_x, train_y = get_data_pickle(path)

    params['name'] = str(p)
    model = Model(params)
    model.plot_data(train_x)
    model.train_model(train_x,train_y)

    start_day = 20
    model.iterative_prediction(train_x,start_day,100)
    # if model == 'GP':
    #     from GP_rrn import train_model
    #     train_model(country,pop,data)

    # elif model == 'RNN':
    #     from rrn import train_model
    #     train_model(country,pop,data)










    