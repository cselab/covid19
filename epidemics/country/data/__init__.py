from pathlib import Path

COUNTRY_DATA_DIR = Path(__file__).parent

NORMALIZED_COUNTRY_NAMES = {
    # https://github.com/samayo/country-json/blob/master/src/country-by-population.json
    "holy see (vatican city state)": "Vatican City",
    "hongkong": "Hong Kong",
    "libryan arab jamahiriya": "Libya",
    "russia": "Russian Federation",

    # hgis
    "uk": "United Kingdom",
}

def _normalize_country_name(name):
    return NORMALIZED_COUNTRY_NAMES.get(name.lower(), name)

def country_to_key(name):
    """Convert country name to a unified format."""
    return _normalize_country_name(name).lower()
