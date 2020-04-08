import math

class Values:
    S = 0
    E = 1
    Ir = 2
    Iu = 3
    N = 4


def flatten(matrix):
    """
    >>> flatten([[10, 20, 30], [40, 50]])
    [10, 20, 30, 40, 50]
    """
    return [
        value
        for row in matrix
        for value in row
    ]


def filter_out_nans_wrt(a, b):
    """
    >>> filter_out_nans_wrt([10, 20, 30, 40], [100, nan, 300, nan])
    ([10, 30], [100, 300])
    """
    assert len(a) == len(b), (len(a), len(b))
    out_a = []
    out_b = []
    for i in range(len(b)):
        if not math.isnan(b[i]):
            out_a.append(a[i])
            out_b.append(b[i])
    return out_a, out_b


def flatten_and_remove_nans(matrix):
    """
    >>> flatten_and_remove_nans([[10, nan, 20], [nan, 30]])
    [10, 20, 30]
    """
    return [
        value
        for row in matrix
        for value in row
        if not math.isnan(value)
    ]
