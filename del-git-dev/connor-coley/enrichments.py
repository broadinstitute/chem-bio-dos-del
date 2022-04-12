from scipy import stats
import numpy as np
from scipy.optimize import fsolve
from functools import partial
from scipy.optimize import fsolve, minimize, minimize_scalar

# FOLLOW: Asymptotic tests based on normal approximations
# https://sci-hub.tw/https://onlinelibrary.wiley.com/doi/pdf/10.1002/bimj.200710403

# Copy code from https://github.com/statsmodels/statsmodels/issues/2718
# copied from statsmodels.stats.weightstats
def _zstat_generic2(value, std_diff, alternative):
    '''generic (normal) z-test to save typing

    can be used as ztest based on summary statistics
    '''
    zstat = value / std_diff
    if alternative in ['two-sided', '2-sided', '2s']:
        pvalue = stats.norm.sf(np.abs(zstat))*2
    elif alternative in ['larger', 'l']:
        pvalue = stats.norm.sf(zstat)
    elif alternative in ['smaller', 's']:
        pvalue = stats.norm.cdf(zstat)
    else:
        raise ValueError('invalid alternative')
        
    pvalue = sys.float_info.min if pvalue == 0 else pvalue
    return zstat, pvalue

def poisson_twosample(count1, exposure1, count2, exposure2, ratio_null=1,
                      method='sqrt', alternative='smaller', ret_stat=False):
    '''test for ratio of two sample Poisson intensities

    If the two Poisson rates are g1 and g2, then the Null hypothesis is

    H0: g1 / g2 = ratio_null

    against one of the following alternatives

    H1_2-sided: g1 / g2 != ratio_null
    H1_larger: g1 / g2 > ratio_null
    H1_smaller: g1 / g2 < ratio_null

    Parameters
    ----------
    count1: int
        Number of events in first sample
    exposure1: float
        Total exposure (time * subjects) in first sample
    count2: int
        Number of events in first sample
    exposure2: float
        Total exposure (time * subjects) in first sample
    ratio: float
        ratio of the two Poisson rates under the Null hypothesis. Default is 1.       
    method: string
        Method for the test statistic and the p-value. Defaults to `'score'`.
        Current Methods are based on Gu et. al 2008
        Implemented are 'wald', 'score' and 'sqrt', see Notes
    alternative : string
        The alternative hypothesis, H1, has to be one of the following

           'two-sided': H1: ratio of rates is not equal to ratio_null (default)
           'larger' :   H1: ratio of rates is larger than ratio_null
           'smaller' :  H1: ratio of rates is smaller than ratio_null

    Returns
    -------
    stat, pvalue two-sided 

    not yet
    #results : Results instance
    #    The resulting test statistics and p-values are available as attributes.


    Notes
    -----
    'wald': method W1A, wald test, variance based on separate estimates
    'score': method W2A, score test, variance based on estimate under Null
    'wald-log': W3A
    'score-log' W4A
    'sqrt': W5A, based on variance stabilizing square root transformation

    References
    ----------
    Gu, Ng, Tang, Schucany 2008: Testing the Ratio of Two Poisson Rates,
    Biometrical Journal 50 (2008) 2, 2008

    '''

    # shortcut names
    y1, n1, y2, n2 = count1, exposure1, count2, exposure2
    d = n2 / n1
    r = ratio_null
    r_d = r / d

    if method in ['score']:
        stats = (y1 - y2 * r_d) / np.sqrt((y1 + y2) * r_d)
    elif method in ['wald']:
        stats = (y1 - y2 * r_d) / np.sqrt(y1 + y2 * r_d**2)
    elif method in ['sqrt']:
        stats = 2 * (np.sqrt(y1 + 3 / 8.) - np.sqrt((y2 + 3 / 8.) * r_d))
        stats /= np.sqrt(1 + r_d)
        
    if ret_stat:
        return stats

    return _zstat_generic2(stats, 1, alternative)

def R_from_z(k2, n2, k1, n1, z):
    '''Analytical solution for getting R from a fixed z value. Uses sqrt method'''
    a = np.power(z, 2) / 4 - (k2 + 3/8)
    b = 2 * np.sqrt(k1 + 3/8) * np.sqrt(k2 + 3/8)
    c = np.power(z, 2) / 4 - (k1 + 3/8)
    x = (-b - np.sign(z) * np.sqrt(np.clip(np.power(b, 2) - 4*a*c, 0, np.inf))) / (2*a)
    return np.power(x, 2) * n2/n1

def R_range(k1, n1, k2, n2, zstat=2):
    return (
        R_from_z(k1, n1, k2, n2, 0),
        R_from_z(k1, n1, k2, n2, -zstat),
        R_from_z(k1, n1, k2, n2, +zstat)
    )

def R(cts, zstat=0):
    k1, n1, k2, n2 = cts
    return R_from_z(k2, n2, k1, n1, zstat)

def R_lb(cts):
    return R(cts, -2)

def R_ub(cts):
    return R(cts, +2)

def R_ranges(k1s, n1s, k2s, n2s, zstat=2, **kwargs):
    return R_range(k1s, n1s, k2s, n2s, zstat)

if __name__ == '__main__':


    c1s = np.array([10, 20, 10, 40]*10000)
    e1s = np.ones(c1s.shape) * c1s.sum()

    c2s = np.array([20, 10, 10, 40]*10000)
    e2s = np.ones(c2s.shape) * c1s.sum()

    (Rs, Rs_lb, Rs_ub) = R_ranges(c1s, e1s, c2s, e2s)
    print(Rs[-1])
    print(Rs_lb[-1])
    print(Rs_ub[-1])




