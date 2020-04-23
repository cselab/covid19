from scipy import interpolate
import numpy as np

def preprocess_data(cases):
    # data is a dictionary of lists, each list the times seris of a field
    # 1. Replaces all leading nans with zeros
    # 2. Truncates end of all time series if missing data in confirmed
    # 3. Interpolates gaps in data (when nan)

    cases.confirmed, len_confirmed       = preprocess_field(np.array(cases.confirmed),'confirmed')
    cases.recovered, len_recovered       = preprocess_field(np.array(cases.recovered),'recovered')
    cases.deaths, len_deaths             = preprocess_field(np.array(cases.deaths),'deaths')

    cases.icu, len_icu                   = preprocess_field(np.array(cases.icu),'icu')
    cases.hospitalized, len_hospitalized = preprocess_field(np.array(cases.hospitalized),'hospitalized')
    cases.ventilated, len_ventilated     = preprocess_field(np.array(cases.ventilated),'ventilated')
    cases.released, len_released         = preprocess_field(np.array(cases.released),'released')

    # For now truncate at number of confirmed
    cases.confirmed = truncate_field(cases.confirmed,len_confirmed)
    cases.recovered = truncate_field(cases.recovered,len_confirmed)
    cases.deaths = truncate_field(cases.deaths,len_confirmed)
    
    cases.icu = truncate_field(cases.icu,len_confirmed)
    cases.hospitalized = truncate_field(cases.hospitalized,len_confirmed)
    cases.ventilated = truncate_field(cases.ventilated,len_confirmed)
    cases.released = truncate_field(cases.released,len_confirmed)

    return cases

def truncate_field(dat,len):
    if dat is not None:
        return dat[:len]
    else:
        return None

def preprocess_field(dat,field):
    print('Processing {}'.format(field))

    if  np.isnan(dat).all():
        print('No data available')
        dat = None
        length = None
    else:
        n = len(dat)
        # 1. Replace all leading nan with zeros
        first_day = np.where(np.isfinite(dat))[0][0]
        dat[:first_day] = 0

        # 2. Remove nan data at the end
        last_day = np.where(np.isfinite(dat))[0][-1]
        dat = dat[:last_day+1]
        print('Removing last {} days in {}'.format(n-last_day,field))

        # 3. Interpolate missing data
        nan_idx = np.where(np.isnan(dat))[0]
        if len(nan_idx) > 0:
            intervals = sublist_split(nan_idx)
            for interval in intervals:
                a = dat[interval[0]-1]
                b = dat[interval[-1]+1]
                dat[interval] = interpolate_missing_data(a,b,len(interval))
        length = len(dat)
        print('{} days of data available for {}'.format(length,field))

    return dat, length

def interpolate_missing_data(a,b,n):
    # for now linear needs to be changed
    x = np.linspace(0,1,n+2)
    f = interpolate.interp1d([0,1], [a,b])
    return f(x[1:-1])

def sublist_split(nan_idx): 
    idx = np.where(np.diff(nan_idx)!=1)[0]
    intervals = []
    previous_interval_length = 0
    for i in idx:
        i -=previous_interval_length
        intervals.append(nan_idx[:i+1])
        previous_interval_length = sum((len(interval) for interval in intervals))
        nan_idx = nan_idx[i+1:]

    intervals.append(nan_idx)
    return intervals