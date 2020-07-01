import numpy as np

from scipy.stats import hypergeom
from scipy.optimize import curve_fit

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
stats = importr('stats')


def false_discovery(pvalues, alpha=0.05):
    """Benjamini-Hochberg procedure for controlling false discovery rate
    """
    pvalues = np.array(pvalues)
    sorter = np.argsort(pvalues)
    n = len(pvalues)
    sig = np.zeros(n).astype(bool)
    for idx, pvalue in enumerate(pvalues[sorter]):
        if pvalue <= alpha * (idx + 1) / n:
            sig[sorter[idx]] = True
        else:
            break

    return sig


def _odds_ratio(table, zero_correction=True):
    """Computes odds ratio from 2x2 contingency table

    [[a, b],
     [c, d]]

    Uses Haldane-Anscombe correction (substitutes 0.5 for 0 values of
    b or c) if zero_correction is set to True.
    """
    ((a, b), (c, d)) = table + zero_correction * 0.5
    se = np.sqrt(np.sum([
        (1/a) + (1/b) + (1/c) + (1/d)
    ]))
    return (a * d) / (b * c), se


def fisher_exact(table, side="two.sided", zero_correction=True):
    """Computes fisher exact odds ratio.
    
    Output is almost exactly the same as scipy.stats.fisher_exact but here allows for
    using Haldane–Anscombe correction (substitutes 0.5 for 0 values in the table, whereas
    the scipy.stats version and R version fisher.test use integers only).
    
    For 95% confidence interval, uses confidence intervals computed by R function fisher.test
    """
    if side not in ("greater", "less", "two.sided"):
        raise ValueError("side parameter must be one of 'greater', 'less', or 'two.sided'")

    # Compute the p value
    # For all possible contingency tables with the observed marginals, compute the hypergeom
    # pmf of that table. Sum the p of all tables with p less than or equal to the hypergeom
    # probability of the observed table.
    N = np.sum(table)
    K = np.sum(table[:, 0])
    n = np.sum(table[0])

    odds_ratio, se = _odds_ratio(table, zero_correction=zero_correction)

    a_min = np.max([0, table[0][0] - table[1][1]])
    a_max = np.min([K, n])
    
    p_observed = hypergeom(N, K, n).pmf(table[0][0])
    p_value = 0.0
    for a in np.arange(a_min, a_max + 1):
        possible_table = np.array([
            [a, n - a],
            [K - a, N - n - K + a]
        ])
        p = hypergeom(N, K, n).pmf(a)
        
        if side == "greater":
            if _odds_ratio(possible_table)[0] >= odds_ratio:
                p_value += p
        elif side == "less":
            if _odds_ratio(possible_table)[0] <= odds_ratio:
                p_value += p
        elif side == "two.sided":
            if p <= p_observed:
                p_value += p

    if side == "greater":
        interval95 = [np.exp(np.log(odds_ratio) - (1.645 * se)), np.inf]
    elif side == "less":
        interval95 = [0, np.exp(np.log(odds_ratio) + (1.645 * se))]
    elif side == "two.sided":
        interval95 = [
                np.exp(np.log(odds_ratio) - (1.96 * se)),
                np.exp(np.log(odds_ratio) + (1.96 * se))
        ]

    return odds_ratio, np.array(interval95), p_value


def r_fisher(table, side="two.sided", zero_correction=True):
    # Get 95% confidence interval from R function fisher.test
    # Use a table with 1 added to zeros if zero_correction is on
    # (this is just for the confidence interval)
    ci_table = table.copy()
    if zero_correction:
        ci_table[ci_table == 0] += 1
    v = robjects.IntVector(np.array(ci_table).flatten())
    m = robjects.r['matrix'](v,nrow=2)
    r_result = stats.fisher_test(m, alternative=side)

    return r_result[2][0], np.array(r_result[1]), r_result[0][0]


def jackknife(samples, estimator):
    """Compute standard error of statistic on given samples

    samples: numpy array of sampled values
    estimator: function that takes numpy array and estimates some statistic (e.g. np.mean)
    
    Returns estimate of standard error of estimator
    """
    jk_n = []
    for i in range(len(samples)):
        jk_n.append(estimator(np.concatenate([samples[:i], samples[i + 1:]])))
    
    jk_n = np.array(jk_n)
    jk_all = estimator(np.array(samples))
    n = len(samples)

    se = np.sqrt(((n - 1) / n) * np.sum((jk_n - jk_all) ** 2))
    return se


def get_odds_ratio_matrix(group1, group2, key):
    """Generate contingency matrix of an in group response and out of group response columns

    |         group1         |         group2         |
    |------------------------|------------------------|
    | #(group1[key] == True) | #(group2[key] == True) |
    | #(group1[key] != True) | #(group2[key] != True) |

    """
    contingency_table = [
        [len(group1[group1[key] == True]),
        len(group2[group2[key] == True])],
        [len(group1[group1[key] == False]),
        len(group2[group2[key] == False])]
    ]

    return np.array(contingency_table)


def compute_odds_ratio(
        group,
        versus,
        zero_correction=True,
        side="two.sided",
    ):
    """Compute odds ratio on an in group and out group
   
    group and versus are pandas DataFrame objects representing
    trials from two conditions. They each should have a boolean column
    named "Response" indicating behavioral response.
    """
    table = get_odds_ratio_matrix(group, versus)    
    odds, interval, pvalue = fisher_exact(table, side=side)

    return odds, interval, pvalue


def linreg(x, y):
    """Perform a simple linear regression on x, y arrays

    Returns:
        popt: optimal values of the parameters (a, b)
        pcov: estimated covariance of the estimated values of popt
        fit_fn: best fit line function, with parameters popt already filled in
        r_squared: R squared value
        r_adj: adjusted R squared value
        sigma_ab: standard deviation of best fit values in popt (squart root of diagonal of cov)
    """
    def lin(x, a, b):
        return x * a + b

    popt, pcov = curve_fit(lin, x, y)
    sigma_ab = np.sqrt(np.diagonal(pcov))
    residuals = y - lin(x, *popt)
    
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    n = len(x)
    k = 1
    r_adj = 1 - ((1 - r_squared) * (n-1) / (n-k-1))

    def fit_fn(x):
        return lin(x, *popt)
    
    return popt, pcov, fit_fn, r_squared, r_adj, sigma_ab

