import numpy as np
from multiprocessing import Pool

class WorkerPool(object):
    def __init__(self, cores):
        self.pool = Pool(processes=cores)
        self.size = cores
    def map(self, function, tasks):
        return self.pool.map(function, tasks)

def priorTransformFromJs(p,js):
    pt = np.zeros(len(p))
    for idx in range(len(p)):
        lb = js['Distributions'][idx]['Minimum']
        ub = js['Distributions'][idx]['Maximum']

        pt[idx] = lb+p[idx]*(ub-lb)
    return pt



def getPosteriorFromResult(result):
    
    from dynesty import utils as dyfunc

    weights = np.exp(result.logwt - result.logz[-1]) # normalized weights
    samples = dyfunc.resample_equal(result.samples, weights)
    return samples
