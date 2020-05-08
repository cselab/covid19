import epidemics.data.swiss_cantons as swiss_cantons
import numpy as np

class Data:
    total_infected = None
    time = None
    population = None
    names = None
    commute_matrix = None  # Cij
    commute_airports = None  # Qa
    commute_borders = None  # Qb

    """
    Fit importance.
    `numpy.ndarray()`, (n_regions,)
    defaults to 1, larger values increase the importance of a region
    """
    fit_importance = None

    """
    Indices of regions for beta correction.
    `dict(varname -> list[region_index])`
    `varname` is ("beta_corr0", ...)
    """
    beta_corr_regions = dict()

def moving_average(x, w):
    """
    x: `numpy.ndarray`, (N)
    w: int
      Window half-width.
    Returns:
    xa: `numpy.ndarray`, (N)
      Array `x` averaged over window [-w,w].
    """
    s = np.zeros_like(x)
    q = np.zeros_like(x)
    for i in range(len(x)):
        for j in range(max(0, i - w), min(i + w + 1, len(x))):
            if np.isfinite(x[j]):
                s[i] += x[j]
                q[i] += 1
    xa = s / q
    return xa


def fill_nans_nearest(x):
    """
    x: `numpy.ndarray`, (N)
    Returns:
    `numpy.ndarray`, (N)
      Array `x` where each NaN value is replaced by a finite value from:
      - nearest index to the left
      - if not found, nearest index to the right
    """
    def left(x, i):
        for v in reversed(x[:i + 1]):
            if np.isfinite(v):
                return v
        return x[i]

    def right(x, i):
        for v in x[i:]:
            if np.isfinite(v):
                return v
        return x[i]

    x = np.array(x)
    for i in range(len(x)):
        x[i] = left(x, i)
        x[i] = right(x, i)
    return x


#
def fill_nans_interp(t, x):
    from scipy.interpolate import interp1d
    x = np.copy(x)
    nans = np.isnan(x)

    f = interp1d(t[~nans], x[~nans], fill_value='extrapolate')
    for i in range(len(x)):
        if np.isnan(x[i]):
            x[i] = f(t[i])
    return x


def get_data_switzerland() -> Data:
    data = Data()
    from epidemics.data.combined import RegionalData
    regionalData = RegionalData('switzerland')
    data.total_infected = regionalData.infected

    data.total_infected = moving_average(data.total_infected, 0)
    data.time = regionalData.time
    data.total_infected = data.total_infected
    data.time = data.time
    data.population = regionalData.populationSize
    return data


def get_all_canton_keys() -> list:
    key_to_population = swiss_cantons.CANTON_POPULATION
    keys = sorted(list(key_to_population))
    return keys


def get_commute_matrix(keys):
    """
    Returns:
    Cij: `numpy.ndarray`, (len(keys), len(keys))
    Cij[work][home] is the number of people registered in canton
    `home` and working in `work`. Data from bfs.admin.ch (2014).
    """
    json = swiss_cantons.get_Cij_home_work_bfs()
    out = np.zeros((len(keys), len(keys)))
    for i1, c1 in enumerate(keys):
        for i2, c2 in enumerate(keys):
            out[i1][i2] = json[c1][c2]
    return out


def get_data_switzerland_cantons(keys) -> Data:
    """
    keys: `array_like` or None
    List of canton keys to select (e.g ['ZH', 'TI'])
    """
    key_to_total_infected = swiss_cantons.fetch_openzh_covid_data()
    key_to_population = swiss_cantons.CANTON_POPULATION

    n_regions = len(keys)

    data = Data()
    if n_regions == 0:
        return data

    nt = len(key_to_total_infected[keys[0]])
    data.time = np.arange(nt)
    data.total_infected = np.empty((n_regions, nt))
    data.population = np.empty((n_regions))

    for i, k in enumerate(keys):
        Itotal = key_to_total_infected[k]
        #Itotal = fill_nans_nearest(Itotal)
        Itotal[0] = 0
        Itotal = fill_nans_interp(data.time, Itotal)
        #Itotal = moving_average(Itotal, 2) # XXX
        data.total_infected[i, :] = Itotal[:]
        data.population[i] = key_to_population[k]
    data.name = keys
    data.commute_matrix = get_commute_matrix(keys)
    data.commute_airports = get_infected_commuters_airports(keys)
    data.commute_borders = get_infected_commuters_borders(keys)

    sel = -1
    sel = 59  # XXX limit data to prevent `free(): invalid next size (fast)`
    #sel = 40
    #sel = 61  # XXX gives `free(): invalid next size (fast)`
    data.time = data.time[:sel]
    data.total_infected = data.total_infected[:, :sel]

    return data


def get_data_synthetic(ode) -> Data:
    data = Data()
    n_regions = 3
    N = np.array([8e6] * n_regions) / (10**np.arange(n_regions))
    I0 = np.array([16] * n_regions)
    t = np.arange(0, 60)
    params = {'N': N}
    params.update(ode.params_fixed)
    S0 = N - I0
    S, I = ode.solve_SI(params, [0, max(t)], [S0, I0], t)
    data.total_infected = N[:, None] - S
    data.time = t
    data.population = N
    return data


def get_infected_commuters_borders(keys):
    FR = 'france'
    GE = 'germany'
    AU = 'austria'
    IT = 'italy'

    cases = {
        AU: 14710,
        FR: 112606,
        GE: 145184,
        IT: 178972,
    }

    population = {
        AU: 8902600,
        FR: 67076000,
        GE: 83149300,
        IT: 60238522,
    }

    travel = [
        ('ZH', GE, 10404.7),
        ('BE', FR, 3514.7),
        ('LU', GE, 613.8),
        ('UR', IT, 46.0),
        ('SZ', AU, 372.9),
        ('OW', IT, 117.6),
        ('NW', IT, 88.2),
        ('GL', AU, 61.4),
        ('ZG', GE, 996.2),
        ('FR', FR, 1027.9),
        ('SO', FR, 2151.6),
        ('BS', GE, 33932.4),
        ('BL', GE, 22318.4),
        ('SH', GE, 4932.3),
        ('AR', AU, 400.4),
        ('AI', AU, 99.7),
        ('SG', AU, 9199.6),
        ('GR', IT, 6998.4),
        ('AG', GE, 13915.3),
        ('TG', GE, 5586.6),
        ('TI', IT, 67878.4),
        ('VD', FR, 32425.2),
        ('VS', FR, 3079.1),
        ('NE', FR, 12944.3),
        ('GE', FR, 87103.8),
        ('JU', FR, 8640.8),
    ]

    travel = {v[0]: (v[1], v[2]) for v in travel}

    r = np.zeros(len(keys))
    for i, key in enumerate(keys):
        country, commuters = travel[key]
        r[i] = cases[country] * commuters / population[country]
    return r


def get_infected_commuters_airports(keys):
    mlnpass_per_year = {
        "ZH": 31.,
        "GE": 17.5,
        "BS": 8.5,
        "BE": 0.13,
        "TI": 0.09,
        "SG": 0.11,
    }
    pass_per_day = {k: v * 1e6 / 365. for k, v in mlnpass_per_year.items()}

    italy_cases = 197675.
    italy_pop = 60238522.
    prob_infected = italy_cases / italy_pop

    r = np.zeros(len(keys))
    for i, key in enumerate(keys):
        if key in pass_per_day:
            r[i] = pass_per_day[key] * prob_infected
    return r

