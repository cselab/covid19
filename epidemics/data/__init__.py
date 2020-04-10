from pathlib import Path

DATA_DIR = Path(__file__).parent
DATA_CACHE_DIR = DATA_DIR / 'cache'
DATA_FILES_DIR = DATA_DIR / 'files'
DATA_DOWNLOADS_DIR = DATA_DIR / 'downloads'

NORMALIZED_COUNTRY_NAMES = {
    # https://github.com/samayo/country-json/blob/master/src/country-by-population.json
    "holy see (vatican city state)": "Vatican City",
    "hongkong": "Hong Kong",
    "libryan arab jamahiriya": "Libya",
    "russia": "Russian Federation",

    # hgis
    "uk": "United Kingdom",
}
