from epidemics.data import DATA_DOWNLOADS_DIR
from epidemics.data.regions import region_to_key
from epidemics.tools.cache import cache
from epidemics.tools.io import download_and_save

from collections import namedtuple
import datetime

class RegionCasesData:
    """Daily data of number of confirmed cases, recovered cases and death for one region."""
    __slots__ = ('start_date', 'confirmed', 'recovered', 'deaths')

    def __init__(self, start_date, confirmed, recovered, deaths):
        self.start_date = start_date
        self.confirmed = confirmed
        self.recovered = recovered
        self.deaths = deaths

    def __repr__(self):
        return "{}(start_date={}, confirmed={}, recovered={}, deaths={})".format(
                self.__class__.__name__, self.start_date, self.confirmed,
                self.recovered, self.deaths)

    def get_confirmed_at_date(self, date):
        idx = (date - self.start_date).days
        if idx < 0 or idx >= len(self.confirmed):
            return 0
        return self.confirmed[idx]

    def get_date_of_first_confirmed(self):
        for day, value in enumerate(self.confirmed):
            if value:
                return self.start_date + datetime.timedelta(days=day)
        raise Exception("Region does not even have confirmed cases.")


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
    start_date = datetime.date.fromisoformat(rows[0].split(',', 1)[0])
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
                start_date=start_date,
                confirmed=[cell[0] for cell in values],
                recovered=[cell[2] for cell in values],
                deaths=[cell[3] for cell in values])

    return out


def get_region_cases(region):
    data = load_and_process_hgis_data()  # This is cached.
    return data[region_to_key(region)]

