from concurrent.futures import ThreadPoolExecutor
from functools import partial

import numpy as np
from scipy.stats import hypergeom, chi2
from scipy.optimize import curve_fit


def false_discovery(pvalues, alpha=0.05):
    """Benjamini-Hochberg procedure for controlling false discovery rate
    """
    pvalues = np.array(pvalues)
    sorter = np.argsort(pvalues)
    n = len(pvalues)
    sig = np.zeros(n).astype(bool)
    thresholds = alpha * (np.arange(n) + 1) / n

    for i, (pvalue, threshold) in enumerate(zip(pvalues[sorter], thresholds)):
        if pvalue <= threshold:
            sig[sorter[:i + 1]] = True

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

    return odds_ratio, np.array(interval95), p_value, se


def jackknife(samples, estimator, parallel=False, **kwargs):
    """Compute standard error of statistic on given samples

    samples: numpy array of sampled values
    estimator: function that takes numpy array and estimates some statistic (e.g. np.mean)

    Returns estimate of standard error of estimator
    """
    jk_n = []
    n = len(samples)

    # Compute the value of estimator over all n samples
    jk_all = estimator(np.array(samples), **kwargs)

    # Compute value of estimator for each combination of n-1 samples
    map_data = [
        np.concatenate([samples[:i], samples[i+1:]])
        for i in range(len(samples))
    ]
    if parallel:
        with ThreadPoolExecutor(4) as pool:
            jk_n = list(pool.map(partial(estimator, **kwargs), map_data))
    else:
        jk_n = [partial(estimator, **kwargs)(s) for s in map_data]
    jk_n = np.array(jk_n)

    # TODO: Estimating psths with the psueodvalues method comes out really bad
    # I don't know why at the moment so just skip that...

    # Compute pseudo values for samples (in n -> inf limit)
    jk_pseudo_values = [(n * jk_all - (n - 1) * jk_n[i]) for i in range(n)]

    est_mean = np.mean(jk_pseudo_values)
    est_var = (1 / n) * np.var(jk_pseudo_values)
    est_sem = np.sqrt(est_var)

    return est_mean, est_sem


def get_odds_ratio_matrix(group1, group2, key):
    """Generate contingency matrix of an in group response and out of group response columns

    |         group1         |         group2         |
    |------------------------|------------------------|
    | #(group1[key] == True) | #(group2[key] == True) |
    | #(group1[key] != True) | #(group2[key] != True) |

    """
    if key is None:
        contingency_table = [
            [len(group1[group1 == True]),
            len(group2[group2 == True])],
            [len(group1[group1 == False]),
            len(group2[group2 == False])]
        ]
    else:
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
        key=None,
        zero_correction=True,
        side="two.sided",
    ):
    """Compute odds ratio on an in group and out group

    group and versus are pandas DataFrame objects representing
    trials from two conditions. They each should have a boolean column
    named "Response" indicating behavioral response.
    """
    table = get_odds_ratio_matrix(group, versus, key=key)
    odds, interval, pvalue = fisher_exact(table, side=side)

    return odds, interval, pvalue


def bootstrap(func, *args, iters=10000):
    """Return bootstrapped standard error for func with 1+ args
    """
    bootstrap_estimates = []
    for _ in range(iters):
        sampled_args = []
        for arg in args:
            sampled_args.append(
                np.random.choice(arg, replace=True, size=len(arg))
            )
        bootstrap_estimates.append(func(*sampled_args))

    return np.std(bootstrap_estimates)


def likelihood_ratio_test(base_model_result, alternate_model_result):
    """Perform a likelihood ratio to compare nested models

    Parameters
    ----------
    base_model_result : statsmodels.base.model.LikelihoodResultsWrapper
        A results wrapper from statsmodels that is the result of fitting a statsmodels model.
        For example, the result of `smf.mixedlm("y ~ x", ...).fit()`. It should have a
        reference to the model and .llf, the log likelihood of the model
    alternate_model_result : statsmodels.base.model.LikelihoodResultsWrapper
        A results wrapper from statsmodels (see base_model). This should be the results of fitting
        an alternate model that adds additional parameters to the base_model.

    Returns
    -------
    p : float
        The p-value (between 0 and 1), representing the probability that one should accept the
        null hypothesis (that the models explain the data equally well).
    x : float
        The chi-squared value, equal to twice the difference between the log-likelihood of the
        two models
    dof : int
        The degrees of freedom for the chi-squared distribution, equal to the difference in the
        number of parameters between the base model and alternate model
    """
    assert alternate_model_result.model.k_params > base_model_result.model.k_params
    dof = alternate_model_result.model.k_params - base_model_result.model.k_params
    x = 2 * (alternate_model_result.llf - base_model_result.llf)
    p = 1 - chi2.cdf(x, dof)

    return p, x, dof


def two_to_one_tail(p: float, t: float, mean: float = 0, side: str = "high"):
    """Convert a p-value from a two tail to one tail test"""
    if side not in ("high", "low"):
        raise ValueError("side must be high or low")
    
    if (side == "high" and t > mean) or (side == "low" and t < mean):
        return p / 2
    else:
        return 1 - (p / 2)
