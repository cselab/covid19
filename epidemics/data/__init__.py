from pathlib import Path

# Deprecated.
DATA_DIR = Path(__file__).parent

DATA_CACHE_DIR = DATA_DIR / 'cache'

# Deprecated.
DATA_FILES_DIR = DATA_DIR / 'files'

DATA_DOWNLOADS_DIR = DATA_DIR / 'downloads'

CANTONS_DATA_DIR = DATA_DIR.parent / 'cantons' / 'data'
