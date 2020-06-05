import numpy as np

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
