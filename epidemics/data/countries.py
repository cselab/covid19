NORMALIZED_COUNTRY_NAMES = {
    # https://github.com/samayo/country-json/blob/master/src/country-by-population.json
    "holy see (vatican city state)": "Vatican City",
    "hongkong": "Hong Kong",
    "libryan arab jamahiriya": "Libya",
    "russia": "Russian Federation",

    # hgis
    "uk": "United Kingdom",
}

def normalize_country_name(name):
    return NORMALIZED_COUNTRY_NAMES.get(name.lower(), name)
