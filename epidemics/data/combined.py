import numpy as np

from epidemics.data.cases import get_region_cases
from epidemics.data.population import get_region_population
from epidemics.data.preprocessor import preprocess_data

class RegionalData:
    """Stores cases and population date for a given region.

    Strips away leading zeros in the number of cases.
    """
    def __init__(self, region,preprocess=False):
        self.region = region
        self.populationSize = get_region_population(region)
        self.preprocess = preprocess

        cases = get_region_cases(region)

        if self.preprocess==True:
            cases = preprocess_data(cases)

        skip = next((i for i, x in enumerate(cases.confirmed) if x), None)
        if skip is None:
            raise ValueError(f"Region `{region}` has no cases.")
        else:
            print('Removing {} leading zeros, {} days of usable data'.format(skip,
                                                    len(cases.confirmed[skip:])))
        # TODO: 'infected' or 'confirmed'?
        self.infected  = add_attribute(cases.confirmed,skip)
        self.recovered = add_attribute(cases.recovered,skip)
        self.deaths    = add_attribute(cases.deaths,skip)
        
        self.hospitalized = add_attribute(cases.hospitalized,skip)
        self.icu = add_attribute(cases.icu,skip)
        self.ventilated = add_attribute(cases.ventilated,skip)
        self.released = add_attribute(cases.released,skip)

        self.time = np.asarray(range(len(self.infected)))

def add_attribute(data,skip):
    if data is not None:
        return np.asarray(data[skip:])
    else:
        return None