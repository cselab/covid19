import numpy as np


class SyntheticData:
    """Reads and stores synthetic data.

    """
    def __init__(self, datafile,preprocess=False):
        print("[Epidemics] Reading Synthetic data.")

        self.preprocess = preprocess
        self.read_synthetic_data(datafile)
        self.time = np.asarray(range(len(self.infected)))

    def read_synthetic_data(self, filename):
        f = open(filename,"r")
        
        self.region = f.readline()
        self.populationSize = int(f.readline())
        
        T = int(f.readline())
        self.infected = np.zeros(T)
        for idx in range(T):
            self.infected[idx] = float(f.readline())

        f.close()


