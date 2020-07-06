import numpy as np

from epidemics.data.cases import get_region_cases
from epidemics.data.population import get_region_population
from epidemics.data.preprocessor import preprocess_data, cut_data_intervention

class RegionalData:
    """Stores cases and population date for a given region.

    Strips away leading zeros in the number of cases.
    """
    def __init__(self, region,preprocess=False,up_to_int=False):
        self.region = region
        self.populationSize = get_region_population(region)
        self.preprocess = preprocess
        self.up_to_int = up_to_int

        cases = get_region_cases(region)

        if self.preprocess==True:
            print('Preprocessing')
            cases = preprocess_data(cases)
        elif self.up_to_int==True:
            cases = cut_data_intervention(cases,region)

        threshold = 2e-6*self.populationSize # 2 per million
        skip = next((i for i, x in enumerate(cases.confirmed) if (x>threshold)), None)
        zeros = next((i for i, x in enumerate(cases.confirmed) if x), None)
        
        if skip is None:
            raise ValueError(f"Region `{region}` has no cases.")
        else:
            print('[Epidemics] Data contain {} leading zeros.'.format(zeros))
            print('[Epidemics] Removing {} entries below threshold (2p. million), {} days of usable data'.format(skip,
                                                    len(cases.confirmed[skip:])))

        self.infected  = add_and_cut_attribute(cases.confirmed,skip) # confirmed
        self.recovered = add_and_cut_attribute(cases.recovered,skip)
        self.deaths    = add_and_cut_attribute(cases.deaths,skip)

        self.hospitalized = add_and_cut_attribute(cases.hospitalized,skip)
        self.icu = add_and_cut_attribute(cases.icu,skip)
        self.ventilated = add_and_cut_attribute(cases.ventilated,skip)
        self.released = add_and_cut_attribute(cases.released,skip)

        self.time = np.asarray(range(len(self.infected)))

def add_and_cut_attribute(data,skip):
    if data is not None:
        return np.asarray(data[skip:])
    else:
        return None
