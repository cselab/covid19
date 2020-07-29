from epidemics.data import DATA_DOWNLOADS_DIR
from epidemics.data.cases import RegionCasesData, RegionalDataBase
from epidemics.country.data import COUNTRY_DATA_DIR, country_to_key
from epidemics.country.data.population import get_country_population
from epidemics.utils.cache import cache
from epidemics.utils.date import date_fromisoformat
from epidemics.utils.io import download_and_save

import os
import pandas as pd
import datetime

def sum_china_province_data(data):
    # list compared with wikipedia, incl. hk, without tibet
    # inner mongolian aut. region == neimenggu
    provinces = ['anhui', 'beijing', 'chongqing', 'fujian', 'gansu', 
            'guangdong', 'guangxi', 'guizhou', 'guizhou', 'hainan', 
            'hebei', 'heilongjiang', 'henan', 'hong kong', 'hubei',
            'hunan', 'neimenggu', 'jiangsu', 'jiangxi', 'jilin', 
            'liaoning', 'macau', 'ningxia', 'qinghai', 'shaanxi', 
            'shandong', 'shanghai', 'shanxi', 'sichuan', 'taiwan', 
            'tianjin', 'xinjiang', 'yunnan', 'zhejiang']

    china = data[provinces[0]]
    for province in provinces[1:]:
        china = china + data[province]

    return china


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
        out[country_to_key(name)] = RegionCasesData(
                start_date=start_date,
                confirmed=[cell[0] for cell in values],
                recovered=[cell[2] for cell in values],
                deaths=[cell[3] for cell in values])

    out['china'] = sum_china_province_data(out)

    return out


def get_country_cases(country):
    data = load_and_process_hgis_data()  # Cached
    return data[country_to_key(country)]


def get_lockdown_date(country, path=COUNTRY_DATA_DIR / 'official_lockdowns.csv'):
    """
    Appends country data by columns:
        official_lockdown: date of official lockdown.
    df: `pandas.DataFrame`
        Output of CollectCountryData()
    path: `str`
        Path to csv generated by `official_lockdowns.py`
    """
    csv = pd.read_csv(path)

    df = pd.DataFrame(data={'country': [''] * len(csv), 'date': [''] * len(csv)})

    for i in range(len(csv)):
        df['country'][i] = csv['country'][i]
        df['date'][i] = csv['date'][i].split('[')[0]

    idx = df.index[df['country'] == country.capitalize()].tolist()
    date = df.iloc[idx]['date'].iloc[0]

    print("    Intervention data for {}: {}".format(country, date))
    return datetime.datetime.strptime(date, '%Y-%m-%d').date()


def cut_data_intervention(cases, country):
    print("[Epidemics][Preprocessing]")
    print("    Removing data pre intervention")

    lockdown_date = get_lockdown_date(country)

    n_days = (lockdown_date - cases.start_date).days + 14
    print(f"    Using {n_days} first days")

    # For now truncate at number of confirmed
    cases.confirmed = truncate_field(cases.confirmed,n_days)
    cases.recovered = truncate_field(cases.recovered,n_days)
    cases.deaths = truncate_field(cases.deaths,n_days)

    cases.icu = truncate_field(cases.icu,n_days)
    cases.hospitalized = truncate_field(cases.hospitalized,n_days)
    cases.ventilated = truncate_field(cases.ventilated,n_days)
    cases.released = truncate_field(cases.released,n_days)

    return cases


class CountryData(RegionalDataBase):
    def __init__(self, country, up_to_int, **kwargs):
        population = get_country_population(country)
        cases  = get_country_cases(country)
        intDay = get_lockdown_date(country)
        if up_to_int:
            cases = cut_data_intervention(cases, region)
        super().__init__(region=country, populationSize=population, cases=cases, \
                interventionDay=intDay, **kwargs)
        self.up_to_int = up_to_int
