import sys
import numpy as np


class SyntheticData:
    """Reads and stores synthetic data.

    """
    def __init__(self, datafile, useInfected, useDeaths, preprocess=False):
        print("[Epidemics] Reading Synthetic data.")

        self.preprocess = preprocess
        self.read_synthetic_data(datafile, useInfected, useDeaths)
        self.time = np.asarray(range(len(self.infected)))

    def read_synthetic_data(self, filename, useInfected, useDeaths):
        f = open(filename,"r")
        
        self.region = f.readline()
        self.populationSize = int(f.readline())
        
        T = int(f.readline())

        if (useInfected == True and useDeaths == False):
            self.infected = np.zeros(T)
            self.deaths   = np.zeros(T)
            self.tact     = None        # dummy
            for idx in range(T):
                self.infected[idx] = float(f.readline())
        
        elif (useInfected == True and useDeaths == True):
            self.infected = np.zeros(int(T/2))
            self.deaths   = np.zeros(int(T/2))
            self.tact     = None        # dummy
            for idx in range(int(T/2)):
                self.infected[idx] = float(f.readline())
            for idx in range(int(T/2)):
                self.deaths[idx]   = float(f.readline())

        else:
            print("choose either (True/True) or (True/False) for (--useInfected/--useDeaths), exit..",flush=True)
            sys.exit()

        f.close()


