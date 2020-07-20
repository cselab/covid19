from epidemics.country.data import COUNTRY_DATA_DIR, country_to_key
from epidemics.utils.cache import cache

import json

def load_samayo_population_data():
    """
    Returns a dictionary {country name: population}.

    Source: https://github.com/samayo/country-json/blob/master/src/country-by-population.json
    """

    with open(COUNTRY_DATA_DIR / 'country-by-population.json') as f:
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
    countries = {country_to_key(name): population for name, population in raw.items()}

    countries[country_to_key('Kosovo')] = 1810463
    countries[country_to_key('Montenegro')] = 631219
    countries[country_to_key('Serbia')] = 6963764
    countries[country_to_key('US')]     = countries['united states']
    return countries


def get_country_population(country):
    """Return the population of the given country."""
    data = get_population_of_all_countries()  # Cached.
    return data[country_to_key(country)]
