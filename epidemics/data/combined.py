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
        self.preprocess =preprocess

        cases = get_region_cases(region)

        if self.preprocess==True:
            cases.data = preprocess_data(cases.data)

        skip = next((i for i, x in enumerate(cases.data['confirmed']) if x), None)
        if skip is None:
            raise ValueError(f"Region `{region}` has no cases.")

        # TODO: 'infected' or 'confirmed'?
        self.infected  = add_attribute(cases.data,'confirmed',skip)
        self.recovered = add_attribute(cases.data,'recovered',skip)
        self.deaths    = add_attribute(cases.data,'deaths',skip)
        
        self.hospitalized = add_attribute(cases.data,'hospitalized',skip)
        self.icu = add_attribute(cases.data,'icu',skip)
        self.ventilated = add_attribute(cases.data,'ventilated',skip)
        self.released = add_attribute(cases.data,'released',skip)

        self.time = np.asarray(range(len(self.infected)))

def add_attribute(data,field,skip):
    if field in data.keys() and not np.isnan(data[field]).all():
        return data[field][skip:]
    else:
        return np.nan