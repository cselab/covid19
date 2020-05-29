#!/usr/bin/env python

import json
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
# from european_countries import EUROPEAN_COUNTRIES
from model import Model
presets = "./countries.json"
import argparse

def filter_few(confirmed, population):
    i = 0
    for i, c in enumerate(confirmed):
        if c > 5 and c > 2e-6 * population:
            break
    return confirmed[i:]

def filter_zeros(confirmed):
    i = 0
    for i, c in enumerate(confirmed):
        if c > 0:
            break
    return confirmed[i:]

def get_data(country):

    with open('./data/countries.json') as f:
        js = json.loads(f.read())

    js_data = [js_i for js_i in js if js_i['country'] in country]
    js_data = js_data[0]

    assert(country == js_data["country"])
    data = js_data["confirmed"]
    pop = js_data["population"]

    data = filter_zeros(data)

    data = np.array(data)
    data = data/pop
    train_input  = data[:-1]
    train_target = data[1:] - data[:-1]

    train_x = np.expand_dims(train_input, axis=1)
    train_y = np.expand_dims(train_target, axis=1)

    return train_x, train_y

if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--country', '-c', default='Switzerland')
    parser.add_argument('--n_features', '-nf', default=1,type=int)
    parser.add_argument('--seq_length', '-sl', default=10,type=int)
    parser.add_argument('--hidden_dim', '-hd', default=20,type=int)
    parser.add_argument('--n_layers', '-nl', default=1,type=int)
    parser.add_argument('--n_epochs', '-ne', default=20000,type=int)
    parser.add_argument('--lr', '-lr', default=0.001,type=int)

    params = vars(parser.parse_args())

    train_x, train_y = get_data(params['country'])

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










    