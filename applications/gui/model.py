# Created by Petr Karnakov on 25.05.2020
# Copyright 2020 ETH Zurich

import os
import numpy as np

from epidemics.epidemics import EpidemicsBase
from epidemics.data.combined import RegionalData
from epidemics.tools.tools import save_file


class EpidemicsCountry(EpidemicsBase):
    def __init__(self, **kwargs, futureDays=0, nSamples=2000, nPropagation=100):
        super().__init__(**kwargs)

        self.__process_data()

    def save_data_path(self):
        return (self.dataFolder, self.country, self.modelName)

    def __process_data(self):
        y = self.regionalData.infected
        t = self.regionalData.time
        N = self.regionalData.populationSize
        I0 = y[0]
        S0 = N - I0
        y0 = S0, I0

        self.data['Model']['x-data'] = t[1:]
        self.data['Model']['y-data'] = np.diff(y[0:])

        self.data['Model']['Initial Condition'] = y0
        self.data['Model'][
            'Population Size'] = self.regionalData.populationSize

        T = np.ceil(t[-1] + self.futureDays)
        self.data['Propagation']['x-data'] = np.linspace(0, T, int(T + 1))

        save_file(self.data, self.saveInfo['inference data'],
                  'Data for Inference', 'pickle')
