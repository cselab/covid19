from scipy import interpolate
import numpy as np

def preprocess_data(data):
    # data is a dictionary of lists, each list the times seris of a field
    # 1. Replaces all leading nans with zeros
    # 2. Truncates end of all time series if missing data in infected
    # 3. Interpolates gaps in data (when nan)
    lengths = {}
    for field in data.keys():
        data[field], lengths[field]= preprocess_field(np.array(data[field]),field)
    print(lengths)
    
    # For now truncate at number of infected
    for field in data.keys(): 
        if not np.isnan(data[field]).all(): 
            data[field] = data[field][:lengths['confirmed']]
    return data

def preprocess_field(dat,field):
    print('Processing {}'.format(field))

    if  np.isnan(dat).all():
        print('No data available')

        dat = np.nan
        length = np.nan
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
        previous_interval_length = len(intervals[-1])
        nan_idx = nan_idx[i+1:]

    intervals.append(nan_idx)
    return intervals