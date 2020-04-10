from epidemics.data import DATA_FILES_DIR
from epidemics.data.countries import normalize_country_name 
from epidemics.tools.cache import cache
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
    countries = {}
    for name, population in raw.items():
        name = normalize_country_name(name)
        key = name.lower()
        countries[key] = population

    countries['kosovo'] = 1810463
    countries['montenegro'] = 631219
    countries['serbia'] = 6963764
    return countries


def get_country_population(country):
    """Return the population of the given country."""
    return get_population_of_all_countries()[country.lower()]
