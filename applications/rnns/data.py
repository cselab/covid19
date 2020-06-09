import json, argparse, pickle
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt

from make_synthetic import sir_int_r0

def get_data(type,**kwargs):

    if type == 'country':
        return get_data_json(country=kwargs['country'])
    
    if type == 'pickle':
        return get_data_pickle(path=kwargs['path'])
    
    if type == 'synthetic':
        return generate_data(p=kwargs['p'],T=kwargs['T'])

def get_data_json(country):

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

def get_data_pickle(path):

    dict_data = pickle.load( open(path, "rb" ) )
    data_all = dict_data['Data']
    params = data_all[0][0]
    fields = data_all[0][1]
    data = data_all[0][2]

    pop = dict_data['Population']

    data = np.array(data)
    data = data/pop
    train_input  = data[:,:-1]
    train_target = data[:,1:] - data[:,:-1]

    train_x = np.expand_dims(train_input, axis=2)
    train_y = np.expand_dims(train_target, axis=2)

    return train_x[1,:], train_y[1,:]

def generate_data(p,T):

    N = 1e7
    I0 = 1
    R0 = 0
    S0 = N - I0 - R0
    y0 = (S0, I0, R0)  
    t_eval = np.linspace(0, T, num=T+1)

    out, fields = sir_int_r0(y0, t_eval, N, p )

    cases = -np.diff(out[0])
    cases[np.where(cases<0)] == 0
    infected = np.cumsum(cases)

    data = infected
    train_input  = data[:-1]
    train_target = data[1:] - data[:-1]

    train_x = np.expand_dims(train_input, axis=1)/np.max(train_target)
    train_y = np.expand_dims(train_target, axis=1)/np.max(train_target)

    return train_x, train_y

# ---------------------------------------------------------------------------- #
# ---------------------------------- Utils ----------------------------------- #
# ---------------------------------------------------------------------------- #

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
