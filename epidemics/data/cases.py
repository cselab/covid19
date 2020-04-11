from epidemics.data import DATA_DOWNLOADS_DIR
from epidemics.data.regions import region_to_key
from epidemics.tools.cache import cache
from epidemics.tools.io import download_and_save

from collections import namedtuple

RegionCasesData = namedtuple('RegionCasesData', ('confirmed', 'recovered', 'deaths'))


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
    data = download_and_save(url, DATA_DOWNLOADS_DIR / 'hgis.virus.csv', cache_duration=3600)
    data = data.decode('utf-8')
    header, *rows = data.split('\n')

    # First column is a date.
    days_cells = [row.split(',')[1:] for row in rows[:-days_to_remove]]
    _, *regions = header.split(',')

    out = {}
    for c, name in enumerate(regions):
        values = []
        for day in range(len(days_cells)):
            cell = days_cells[day][c]
            if cell:
                cell = tuple(int(x or 0) for x in cell.split('-'))
            if not cell or len(cell) != 4:
                cell = (0, 0, 0, 0)
            values.append(cell)

        # The meaning of cell[2] is unknown (it seems to be 0 everywhere).
        out[region_to_key(name)] = RegionCasesData(
                confirmed=[cell[0] for cell in values],
                recovered=[cell[2] for cell in values],
                deaths=[cell[3] for cell in values])

    return out


def get_region_cases(region):
    data = load_and_process_hgis_data()  # This is cached.
    return data[region_to_key(region)]

