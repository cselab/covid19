#!/usr/bin/env python3

"""Solve the ODE and plot the results."""

import numpy as np
import matplotlib.pyplot as plt

import argparse
import os
import sys
import collections

sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', '..'))

import epidemics.data
from epidemics.data.swiss_cantons import CANTON_KEYS_ALPHABETICAL, CANTON_POPULATION
from epidemics.cantons.py.model import \
        get_canton_model_data, get_canton_reference_data, \
        get_municipality_model_data, ModelData
from epidemics.cantons.py.plot import Renderer

import libepidemics

class Level:
    canton = "canton"
    municipality = "municipality"

def get_foreign_commuters():
    FR = 'france'
    GE = 'germany'
    AU = 'austria'
    IT = 'italy'

    cases = {
            AU: 14710,
            FR: 112606,
            GE: 145184,
            IT: 178972,
            }

    travel = [
        ('ZH', GE, 10404.7),
        ('BE', FR, 3514.7),
        ('LU', GE, 613.8),
        ('UR', IT, 46.0),
        ('SZ', AU, 372.9),
        ('OW', IT, 117.6),
        ('NW', IT, 88.2),
        ('GL', AU, 61.4),
        ('ZG', GE, 996.2),
        ('FR', FR, 1027.9),
        ('SO', FR, 2151.6),
        ('BS', GE, 33932.4),
        ('BL', GE, 22318.4),
        ('SH', GE, 4932.3),
        ('AR', AU, 400.4),
        ('AI', AU, 99.7),
        ('SG', AU, 9199.6),
        ('GR', IT, 6998.4),
        ('AG', GE, 13915.3),
        ('TG', GE, 5586.6),
        ('TI', IT, 67878.4),
        ('VD', FR, 32425.2),
        ('VS', FR, 3079.1),
        ('NE', FR, 12944.3),
        ('GE', FR, 87103.8),
        ('JU', FR, 8640.8),
    ]

    r = {}
    for v in travel:
        r[v[0]] = cases[v[1]] * v[2]
    return r


def example_run_seiin(data: ModelData, num_days: int, level):
    """Runs the SEIIN model for some set of parameters and some initial conditions."""

    # Parameters.
    # Li2020, Table 1
    params = libepidemics.cantons.seiin.Parameters(
            beta=1.12, mu=0., alpha=1., Z=3.69, D=3.47, theta=1.36)

    # Initial state.
    N0 = list(data.region_population)
    E0 = [0] * data.num_regions
    IR0 = [0] * data.num_regions
    IU0 = [0] * data.num_regions

    src = np.zeros(data.num_regions)
    if level == Level.canton:
        IR0[data.key_to_index['TI']] = 0  # Ticino.
        IU0[data.key_to_index['TI']] = 0

        # sources from airports
        if data.airports:
            k_air = 1.
            for a in data.airports:
                src[data.key_to_index[a[0]]] = a[2] * k_air

        # sources from foreign commuters
        if data.foreign:
            foreign = get_foreign_commuters()
            print(foreign)
            k_foreign = 1e-8
            for k in data.key_to_index:
                src[data.key_to_index[k]] = foreign[k] * k_foreign


    elif level == Level.municipality:
        IR0[data.key_to_index['MUN-5192']] = 0  # Lugano.
        IU0[data.key_to_index['MUN-5192']] = 0

        # sources from airports
        if data.airports:
            k_air = 1.
            for a in data.airports:
                src[data.key_to_index[a[1]]] = a[2] * k_air

        # sources from foreign commuters
    else:
        assert False
    data.ext_com_Iu = [src]

    S0 = [N - E - IR - IU for N, E, IR, IU in zip(N0, E0, IR0, IU0)]
    y0 = S0 + E0 + IR0 + IU0 + N0

    # Run the ODE solver.
    solver = libepidemics.cantons.seiin.Solver(data.to_cpp())
    y0 = libepidemics.cantons.seiin.State(y0)
    return solver.solve(params, y0, t_eval=range(1, num_days + 1))


def example_run_seii_c(data, num_days):
    """Runs the SEII_C model for some set of parameters and some initial conditions."""
    # Parameters. Note that nu is relative to beta.
    params = libepidemics.cantons.seii_c.Parameters(
            beta=3.0, nu=1.0, alpha=0.6, Z=3.69, D=3.47)

    # Initial state.
    E0 = [0] * data.num_regions
    IR0 = [0] * data.num_regions
    IU0 = [0] * data.num_regions

    # if 'TI' in data.key_to_index:
    #     IU0[data.key_to_index['TI']] = 10  # Ticino.
    # else:
    #     IU0[data.key_to_index['MUN-5192']] = 10  # Lugano.

    S0 = [N - E - IR - IU for N, E, IR, IU in zip(data.region_population, E0, IR0, IU0)]
    y0 = S0 + E0 + IR0 + IU0

    # Run the ODE solver.
    solver = libepidemics.cantons.seii_c.Solver(data.to_cpp())
    y0 = libepidemics.cantons.seii_c.State(y0)
    return solver.solve(params, y0, t_eval=range(1, num_days + 1))


def plot_ode_results_canton(data: ModelData, results, air=True, crossborder=True):
    """Plot results from the ODE.

    Arguments:
        results: A list of State objects.
    """

    def frame_callback(rend):
        Ir_max = np.max([state.Ir() for state in results])
        t = rend.get_frame() * (len(results) - 1) // rend.get_max_frame()

        state = results[t]
        values = {}
        texts = {}
        for i, key in enumerate(data.region_keys):
            N = state.N(data.key_to_index[key])
            S = state.S(data.key_to_index[key])
            Ir = state.Ir(data.key_to_index[key])
            Iu = state.Iu(data.key_to_index[key])
            print("{:02d} {} Ir={:4.1f} Iu={:4.1f}".format(i, key, Ir, Iu))
            #values[key] = Ir / Ir_max * 2
            values[key] = (N - S) * 1e-3
            texts[key] = str(int(Ir))
        rend.set_values(values)
        rend.set_texts(texts)

    airports = None
    if data.airports:
        airports = [v[0] for v in data.airports]

    rend = Renderer(frame_callback, data=data, draw_Mij=True, draw_Cij=False,
            airports=airports)
    rend.save_image()
    rend.save_movie(frames=len(results))

def plot_ode_results_munic(data: ModelData, results, air=True, crossborder=True):
    """Plot results from the ODE.

    Arguments:
        results: A list of State objects.
    """

    from epidemics.data.swiss_municipalities import get_cantons, get_name_and_population
    from pandas import Series

    df = get_cantons()
    key_to_canton = Series(df.canton.values,index=df.key).to_dict()

    df = get_name_and_population()
    key_to_name = Series(df.name.values,index=df.key).to_dict()

    def frame_callback(rend):
        t = rend.get_frame() * (len(results) - 1) // rend.get_max_frame()

        plot_NmS = True

        state = results[t]
        Ir_cantons = collections.defaultdict(float)
        NmS_cantons = collections.defaultdict(float)
        for i, key in enumerate(data.region_keys):
            if key in key_to_canton:
                Ir_cantons[key_to_canton[key]] += state.Ir(data.key_to_index[key])
                NmS_cantons[key_to_canton[key]] += \
                        state.N(data.key_to_index[key]) - state.S(data.key_to_index[key])

        n_matches = 0
        # Municipalities
        if rend.draw_zones:
            name_to_value = {}
            for i, key in enumerate(data.region_keys):
                if key in key_to_canton:
                    j = data.key_to_index[key]
                    if plot_NmS:
                        name_to_value[key_to_name[key]] = state.N(j) - state.S(j)
                    else:
                        name_to_value[key_to_name[key]] = state.Ir(j)
            nn = rend.get_zone_names()
            zone_values = np.zeros(len(nn))
            for i,n in enumerate(nn):
                if n in name_to_value:
                    zone_values[i] = np.clip(name_to_value[n] * 0.1, 0., 1.)
                    n_matches += 1
            rend.set_zone_values(zone_values)
            print("found {:} matches between zone and municipality names".format(n_matches))

        Ir_max = np.max(list(Ir_cantons.values()))
        NmS_max = np.max(list(NmS_cantons.values()))
        texts = {}
        values = {}
        # Cantons
        for i,code in enumerate(Ir_cantons):
            print("{:02d} {} Ir={:4.1f} N-S={:4.1f}".format(i, code, Ir_cantons[code], NmS_cantons[code]))
            if plot_NmS:
                #values[code] = NmS_cantons[code] / NmS_max * 2
                values[code] = NmS_cantons[code] * 1e-3
                texts[code] = str(round(NmS_cantons[code]))
            else:
                #values[code] = Ir_cantons[code] / Ir_max * 2
                values[code] = Ir_cantons[code] * 1e-3
                texts[code] = str(round(Ir_cantons[code]))

        rend.set_values(values)
        rend.set_texts(texts)

    airports = None
    if data.airports:
        airports = [v[1] for v in data.airports]

    rend = Renderer(frame_callback, data=data, draw_Mij=False, draw_Cij=False, draw_zones=True, airports=airports)
    rend.save_image()
    rend.save_movie(frames=len(results))


def plot_timeseries(data: ModelData, results, keys, var='S', ref_data=None):
    key_to_population = dict(zip(data.region_keys, data.region_population))

    VAR_NmS = "N-S"
    VAR_S0mS = "S0-S"
    if not isinstance(keys, list):
        keys = [keys]
    fig = plt.figure(figsize=(6, 4))
    ax = fig.gca()
    ax.set_title("   ".join(
        ["$N_{{{:}}}$={:.0f}K".format(key, key_to_population[key]*1e-3) for key in keys]))
    vardesc = {
            'Ir' : "Infected Reported",
            'Iu' : "Infected Unreported",
            VAR_NmS : "N - S(t)",
            VAR_S0mS : "S(0) - S(t)",
            }
    ax.set_xlabel("days")
    ax.set_ylabel(vardesc[var])
    for key in keys:
        if var == VAR_NmS:
            u = np.array([state.S(data.key_to_index[key]) for state in results])
            u = key_to_population[key] - u
        elif var == VAR_S0mS:
            u = np.array([state.S(data.key_to_index[key]) for state in results])
            u = u[0] - u
        else:
            u = [getattr(state, var)(data.key_to_index[key]) for state in results]
        line, = ax.plot(u,label=key)
        if ref_data:
            u = ref_data.cases_per_country[key]
            ax.plot(u, c=line.get_color(), ls='', marker='.')

    if ref_data:
        ax.set_ylim([1, 10 * np.nanmax(ref_data.cases_per_country[key])])
    ax.set_yscale('log')
    ax.legend()
    varname = {
            'Ir' : "Ir",
            'Iu' : "Iu",
            VAR_NmS : "NmS",
            VAR_S0mS : "S0mS",
            }
    fig.tight_layout()
    fig.savefig("{:}_{:}.pdf".format(varname[var], "_".join(keys)))

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('type', type=str, choices=('video', 'timeseries'), help="Plot type.")
    parser.add_argument('days', type=int, default=50, help="Number of days to evaluate.")
    parser.add_argument('--no-foreign', action='store_true', help="Disable foreign commuters from the model.")
    parser.add_argument('--level', type=str, choices=(Level.canton, Level.municipality), default='canton', help="Level of details.")
    parser.add_argument('--model', type=str, choices=('seiin', 'seii_c'), default='seiin', help="Model.")
    args = parser.parse_args(argv)

    if args.level == Level.canton:
        model_data = get_canton_model_data(include_foreign=not args.no_foreign)
        ref_data = get_canton_reference_data()
    elif args.level == Level.municipality:
        model_data = get_municipality_model_data()
        ref_data = None
    else:
        assert False

    # airport passengers (millions per year)
    model_data.airports = [
            ["ZH",  "MUN-0062", 31.,  ],
            ["GE",  "MUN-6623", 17.5, ],
            ["BS",  "MUN-2701", 8.5,  ],
            ["BE",  "MUN-0861", 0.13, ],
            ["TI",  "MUN-5192", 0.09, ],
            ["SG",  "MUN-3237", 0.11, ],
            ]
    model_data.foreign = True

    #model_data.airports = None
    #model_data.foreign = None

    if args.model == 'seiin':
        results = example_run_seiin(model_data, args.days, args.level)
    else:
        model_data.Mij *= 0.0
        results = example_run_seii_c(model_data, args.days)

    if args.type == 'video':
        if args.level == Level.canton:
            plot_ode_results_canton(model_data, results)
        elif args.level == Level.municipality:
            plot_ode_results_munic(model_data, results)
    elif args.type == 'timeseries':
        keys = ['TI', 'ZH', 'AG']
        #plot_timeseries(model_data, results, keys, var='Ir', ref_data=ref_data)
        plot_timeseries(model_data, results, keys, var="N-S", ref_data=ref_data)
        #plot_timeseries(model_data, results, keys, var="S0-S", ref_data=ref_data)


if __name__ == '__main__':
    main(sys.argv[1:])
