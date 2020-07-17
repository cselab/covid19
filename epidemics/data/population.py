from epidemics.data import DATA_FILES_DIR
from epidemics.data.regions import region_to_key
from epidemics.utils.cache import cache
from epidemics.cantons.data.canton_population import CANTON_POPULATION

import json

def load_samayo_population_data():
    """
    Returns a dictionary {country name: population}.

    Source: https://github.com/samayo/country-json/blob/master/src/country-by-population.json
    """

    with open(DATA_FILES_DIR / 'country-by-population.json') as f:
        countries = json.load(f)
    countries = {
        item["country"]: int(item["population"])
        for item in countries
        if item["population"]
    }
    return countries

@cache
def get_population_of_all_countries():
    """Returns a dictionary {country key: population}."""
    raw = load_samayo_population_data()
    countries = {region_to_key(name): population for name, population in raw.items()}

    countries[region_to_key('Kosovo')] = 1810463
    countries[region_to_key('Montenegro')] = 631219
    countries[region_to_key('Serbia')] = 6963764
    return countries


def get_region_population(region):
    """Return the population of the given region."""
    # Currently we only have country population.
    if region in CANTON_POPULATION.keys():
        return CANTON_POPULATION[region]
    else:
        data = get_population_of_all_countries()
        return data[region_to_key(region)]
