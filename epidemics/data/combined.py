import numpy as np

from epidemics.data.cases import get_region_cases
from epidemics.data.population import get_region_population


class RegionalData:
    """Stores cases and population date for a given region.

    Strips away leading zeros in the number of cases.
    """
    def __init__(self, region):
        self.region = region
        self.populationSize = get_region_population(region)

        cases = get_region_cases(region)

        skip = next((i for i, x in enumerate(cases.confirmed) if x), None)
        if skip is None:
            raise ValueError(f"Region `{region}` has no cases.")

        # TODO: 'infected' or 'confirmed'?
        self.infected  = np.asarray(cases.confirmed[skip:])
        self.recovered = np.asarray(cases.recovered[skip:])
        self.deaths    = np.asarray(cases.deaths[skip:])
        self.time      = np.asarray(range(len(cases.confirmed) - skip))
