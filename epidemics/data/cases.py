from epidemics.data.preprocessor import preprocess_data
from epidemics.utils.date import date_fromisoformat
from epidemics.utils.io import download_and_save

from itertools import zip_longest
from collections import namedtuple
from operator import add
import datetime
import numpy as np


# FIXME: The difference between RegionCasesData and RegionalDataBase (i.e.
#        CountryData and CantonsData) is historical. RegionalCasesData was
#        supposed to store only epidemics-related data, while RegionalData
#        contained combined epidemics data with population and other
#        potentially useful data. In principle, RegionCasesData could be
#        removed.

class RegionCasesData:
    """Daily data of number of confirmed cases, recovered cases and death for one region."""

    def __init__(self, start_date, confirmed, recovered, deaths,
                    hospitalized=None, icu=None, released=None, ventilated=None):
        self.start_date = start_date

        self.confirmed = confirmed
        self.recovered = recovered
        self.deaths    = deaths

        self.hospitalized = hospitalized
        self.icu          = icu
        self.released     = released
        self.ventilated   = ventilated

    def __repr__(self):
        return "{}(start_date={}, confirmed={}, recovered={}, deaths={})".format(
                self.__class__.__name__, self.start_date, self.confirmed,
                self.recovered, self.deaths)

    def __add__(self, o):
        start = min(self.start_date, o.start_date)
        confirmed = [x+y for x,y in zip_longest(self.confirmed, o.confirmed, fillvalue=0)]
        recovered = [x+y for x,y in zip_longest(self.recovered, o.recovered, fillvalue=0)]
        deaths    = [x+y for x,y in zip_longest(self.deaths, o.deaths, fillvalue=0)]
        
        return RegionCasesData(start, confirmed, recovered, deaths)

    def get_confirmed_at_date(self, date):
        idx = (date - self.start_date).days
        if idx < 0 or idx >= len(self.confirmed):
            return 0
        return self.confirmed[idx]

    def get_date_of_first_confirmed(self):
        for day, value in enumerate(self.confirmed):
            if value:
                return self.start_date + datetime.timedelta(days=day)
        raise Exception("Region does not even have confirmed cases.")


class RegionalDataBase:  # Base class, not database.
    """Stores cases and population date for a given region.

    Strips away leading zeros in the number of cases.
    """
    def __init__(self, region, *, populationSize, cases, interventionDay, lastDay=datetime.date.today(), preprocess=False):
        self.region = region
        self.populationSize = populationSize
        self.preprocess = preprocess

        if self.preprocess:
            print('Preprocessing')
            cases = preprocess_data(cases)

        fraction  = 1e-8 # 1 ppm
        threshold = fraction*self.populationSize 
        skip  = next((i for i, x in enumerate(cases.confirmed) if (x>threshold)), None)
        zeros = next((i for i, x in enumerate(cases.confirmed) if x), None)
        end   = min((lastDay - cases.start_date).days, len(cases.confirmed))

        if skip is None:
            raise ValueError(f"Region `{region}` has no cases.")
        else:
            print('[Epidemics] Data contain {} leading zeros.'.format(zeros))
            print('[Epidemics] Removing {} entries below threshold ({}/{} pct.), {} days of usable data'.format(skip, threshold, fraction*100, len(cases.confirmed[skip:])))
            print('[Epidemics] Omitting last {} entries'.format(len(cases.confirmed)-end))

        def cut(array):
            if array is not None:
                return np.asarray(array[skip:end])
            else:
                return None

        self.infected  = cut(cases.confirmed) # confirmed
        self.recovered = cut(cases.recovered)
        self.deaths    = cut(cases.deaths)
        print("Infected")
        print(self.infected)
        print("Deaths")
        print(self.deaths)
        self.hospitalized = cut(cases.hospitalized)
        self.icu = cut(cases.icu)
        self.ventilated = cut(cases.ventilated)
        self.released = cut(cases.released)
        
        self.tact = (interventionDay - cases.start_date).days - skip
        self.time = np.asarray(range(len(self.infected)))
