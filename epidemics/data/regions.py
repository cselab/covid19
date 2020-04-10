NORMALIZED_REGION_NAMES = {
    # https://github.com/samayo/country-json/blob/master/src/country-by-population.json
    "holy see (vatican city state)": "Vatican City",
    "hongkong": "Hong Kong",
    "libryan arab jamahiriya": "Libya",
    "russia": "Russian Federation",

    # hgis
    "uk": "United Kingdom",
}

def normalize_region_name(name):
    return NORMALIZED_REGION_NAMES.get(name.lower(), name)

def region_to_key(name):
    return normalize_region_name(name).lower()
