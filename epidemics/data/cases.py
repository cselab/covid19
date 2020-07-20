from epidemics.data import DATA_DOWNLOADS_DIR
from epidemics.cantons.data.canton_population import CANTON_POPULATION

from epidemics.data.regions import region_to_key
from epidemics.utils.cache import cache
from epidemics.utils.date import date_fromisoformat
from epidemics.utils.io import download_and_save

from itertools import zip_longest
from collections import namedtuple
import datetime
from operator import add

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

@cache
def load_and_process_hgis_data(*, days_to_remove=1):
    """Load and reorganize HGIS data.

    Source:
        https://hgis.uw.edu/virus/assets/virus.csv

    Returns:
        A dictionary {region key: RegionCasesData}
    """
    # Note: if you need US and China data, check the coronavirus GUI repo.

    url = 'https://hgis.uw.edu/virus/assets/virus.csv'
    data = download_and_save(url, DATA_DOWNLOADS_DIR / 'hgis.virus.csv', cache_duration=7200)
    data = data.decode('utf-8')
    header, *rows = data.split('\n')

    # First column is a date.
    start_date = date_fromisoformat(rows[0].split(',', 1)[0])
    days_cells = [row.split(',')[1:] for row in rows[:-days_to_remove]]
    _, *regions = header.split(',')

    out = {}
    for c, name in enumerate(regions):
        values = []
        for day in range(len(days_cells)):
            cell = days_cells[day][c]
            if cell:
                cell = tuple(int(0 if x == 'No data' else x or 0) for x in cell.split('-'))
            if not cell or len(cell) != 4:
                cell = (0, 0, 0, 0)
            values.append(cell)

        # The meaning of cell[2] is unknown (it seems to be 0 everywhere).
        out[region_to_key(name)] = RegionCasesData(
                start_date=start_date,
                confirmed=[cell[0] for cell in values],
                recovered=[cell[2] for cell in values],
                deaths=[cell[3] for cell in values])

    return out

@cache
def get_data_of_all_cantons():
    '''
        Gets confirmed(infected), recovered and deaths from the openzh database
    '''
    fields = ['cases','fatalities','released','hospitalized','icu','vent']

    # Get data for all fields
    data = {}
    for field in fields:
        data[field] = get_field_data_all_cantons(field)

    cantons = list(data['cases'].keys())
    cantons.remove('date')
    # Rearrange data per canton
    out = {}
    for canton in cantons:
        recovered = list(map(add, data['fatalities'][canton], data['released'][canton]))
        out[canton] = RegionCasesData(  start_date=data['cases']['date'],
                                        confirmed=data['cases'][canton],
                                        recovered=recovered,
                                        deaths=data['fatalities'][canton],
                                        hospitalized = data['hospitalized'][canton],
                                        icu = data['icu'][canton],
                                        ventilated=data['vent'][canton],
                                        released = data['released'][canton])
    return out

@cache
def get_field_data_all_cantons(field,cache_duration=1e9):

    url = 'https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_'+field+'_switzerland_openzh.csv'
    file_path = 'covid19_'+field+'_switzerland_openzh.csv'
    path = DATA_DOWNLOADS_DIR / file_path

    raw = download_and_save(url, path, cache_duration=cache_duration)
    rows = raw.decode('utf8').split()
    cantons = rows[0].split(',')[1:-1]  # Skip the "Date" and "CH" cell.

    data = {canton: [] for canton in cantons}
    date = []
    for day in rows[1:]:  # Skip the header.
        cells = day.split(',')[1:-1]  # Skip "Date" and "CH".
        date.append(date_fromisoformat(day.split(',')[0]))
        assert len(cells) == len(cantons), (len(cells), len(cantons))
        for canton, cell in zip(cantons, cells):
            data[canton].append(float(cell or 'nan'))
    data['date'] = date
    return data


def collect_province_data(data):
    # list compared with wikipedia, incl. hk, without tibet
    # inner mongolian aut. region == neimenggu
    provinces = ['anhui', 'beijing', 'chongqing', 'fujian', 'gansu', 
            'guangdong', 'guangxi', 'guizhou', 'guizhou', 'hainan', 
            'hebei', 'heilongjiang', 'henan', 'hong kong', 'hubei',
            'hunan', 'neimenggu', 'jiangsu', 'jiangxi', 'jilin', 
            'liaoning', 'macau', 'ningxia', 'qinghai', 'shaanxi', 
            'shandong', 'shanghai', 'shanxi', 'sichuan', 'taiwan', 
            'tianjin', 'xinjiang', 'yunnan', 'zhejiang']

    print(data.keys(),flush=True)
    china = data[provinces[0]]
    for province in provinces[1:]:
        china = china + data[province]

    data['china'] = china
    return data

def get_region_cases(region):
    if region in CANTON_POPULATION.keys():
        data = get_data_of_all_cantons()
        return data[region]
    else:
        data = load_and_process_hgis_data()  # This is cached.
        if region == "china":
            data = collect_province_data(data)
        return data[region_to_key(region)]

