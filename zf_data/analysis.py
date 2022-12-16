import logging
from typing import Iterable, Union
try:
    from functools import cache, cached_property
except ImportError:
    from cached_property import cached_property
    from functools import lru_cache
    cache = lru_cache(maxsize=None)

import numpy as np
import pandas as pd
import scipy.stats

from .stats import jackknife, fisher_exact


logger = logging.getLogger(__name__)


def count_relative_informative_trials(df: pd.DataFrame) -> pd.DataFrame:
    """Compute a relative informative trials seen column from trial dataframe

    Counts informative trials relative to the start of the given dataframe

    Parameters
    ----------
    df : pd.DataFrame
        A dataframe of trials requiring at least the following columns: "Subject",
        "StimulusVocalizerId", "AbsInformativeTrialsSeen"

    Returns
    -------
    df : pd.DataFrame
        A copy of the input dataframe with an additional column, "RelInformativeTrialsSeen".
        This column contains the informative trial counts are counted relative to the
        first appearance of each (subject, vocalizer_id) pair in the dataframe
    """
    df = df.copy()
    df["RelInformativeTrialsSeen"] = df["AbsInformativeTrialsSeen"]
    for _, subdf in df.groupby(["Subject", "StimulusVocalizerId"]):
        df.loc[
            subdf.index,
            "RelInformativeTrialsSeen"
        ] -= subdf.iloc[0]["AbsInformativeTrialsSeen"]

    return df


class Tsvk:
    """Helper class to compute interruption probabilities and odds ratios

    Initialize with a subset of trials

    Example
    -------
    >>> from zf_data import Tsvk, load
    >>> df = load("TrialData.csv")
    >>> T = Tsvk(
            df[
                (df["VocalizerSet"] == "S1") &
                (df["LesionStage"] == "postlesion") &
                (df["SubjectTreatment"] == "NCM")
            ],
            k_max=10
        )

    # Get avg probability of interruption for a given subject across rewarded vocalizers
    >>> T.re.p(subject="BluBlu3434M", k=0)
    0.234

    # Get avg probability of interruption for each subject across nonrewarded vocalizers
    >>> T.nore.p(k=0)
    array([
        ["BluBlu3434M", 0.234],
        ...
        ["RedRed1212F", 0.123],
    ])

    # Get probabilities of interruption over informative trial bins
    >>> T.re.p()
    array([
        [0, 0.12, 0.05],
        ...
        [10, 0.02, 0.01]
    ])
    
    # Compute the logOR as a function of informative trial with stats
    >>> T.logOR()
    array([
        [0, -1.2, 5.4, 0.564, 0.12],
        ...
        [10, 23.2, 4.4, 0.002, 5.65],
    ])

    Arguments
    ---------
    df : pd.DataFrame
    k_max : int
        Maximium informative trial bin to compute out to.
        Leave blank (None) to include all trials
    """
    def __init__(self, df: pd.DataFrame, k_max: int = None):
        self.df = count_relative_informative_trials(df)
        if k_max is None:
            self.k_max = np.max(self.df["RelInformativeTrialsSeen"])
        else:
            self.k_max = min(k_max, np.max(self.df["RelInformativeTrialsSeen"]))

        if k_max != self.k_max:
            logger.debug(f"Provided trials do not reach k_max={k_max} informative trials; using {self.k_max} instead")
        # self.df = self.df[self.df["RelInformativeTrialsSeen"] <= self.k_max]

        self.subjects = self.df.Subject.unique()
        self.k = np.arange(self.k_max + 1)

        # Create a mapping from each subject to all vocalizers with at least one trial
        self.subject_to_vocalizers = dict(
            (subject, self.df[self.df["Subject"] == subject]["StimulusVocalizerId"].unique())
            for subject in self.subjects
        )

        if not len(self.df):
            raise ValueError("Dataframe passed to Tsvk is empty")

        # Check if dataframe has any rewarded or nonrewarded
        self._includes_both_reward_classes = set(["Rewarded", "Nonrewarded"]) == set(self.df["StimulusClass"].unique())

    def filter(self, condition: Iterable) -> 'Tsvk':
        return Tsvk(self.df[condition], k_max=self.k_max)

    @cached_property
    def re(self) -> 'Tsvk':
        if not self._includes_both_reward_classes:
            raise ValueError("Cannot get .re from Tsvk that doesn't have both Rewarded and Nonrewarded trials. It is probably already filtered.")
        return self.filter(self.df["StimulusClass"] == "Rewarded")

    @cached_property
    def nore(self) -> 'Tsvk':
        if not self._includes_both_reward_classes:
            raise ValueError("Cannot get .nore from Tsvk that doesn't have both Rewarded and Nonrewarded trials. It is probably already filtered.")
        return self.filter(self.df["StimulusClass"] == "Nonrewarded")

    @cache
    def _t_svk(self, subject: str, vocalizer: str, k: Union[int, Iterable[int]]):
        """The number of informative trials in the kth bin for subject and vocalizer

        Definition of $$T^{sv}_{k}$$ in the paper.

        Returns
        -------
        T_svk : int
            Number of trials between the kth (exclusive) and the (k+1)th (inclusive)
            informative trial of vocalizer by subject
        n_bins : int
            Number of valid values for k (e.g. if a subject only had 3 informative trials
            and the last k bin is 3, then there were n_bins=4, even if k=range(0, 10)
            was passed into this function.)
        """
        if np.issubdtype(type(k), np.integer):
            selected_trials = self.df[
                (self.df["Subject"] == subject) &
                (self.df["StimulusVocalizerId"] == vocalizer) &
                (self.df["RelInformativeTrialsSeen"] == k)
            ]
        else:
            selected_trials = self.df[
                (self.df["Subject"] == subject) &
                (self.df["StimulusVocalizerId"] == vocalizer) &
                self.df["RelInformativeTrialsSeen"].isin(k)
            ]

        # Generally, we assume that the last bin ends with a non-interruption
        # interrupts: X X X O X X X O X X O
        # k bin:      0 0 0 0 1 1 1 1 2 2 2
        #
        # and so we estimate the probability of interruption in each bin as (T-1)/T,
        # where T is the length of the bin (e.g. 3/4, 3/4, and 2/3). But if are beyond
        # the last non-interrupted trial, the last bin is open ended...
        # interrupts: X X X O X X X O X X O X X
        # k bin:      0 0 0 0 1 1 1 1 2 2 2 3 3
        #
        # ...and it does not make sense to estimate the last probability in bin 3 as 2/2 since
        # we reached the end of the data. Instead, we estimate it as ((2T+1)-1)/(2T+1), or 4/5
        #
        # This next condition handles this case by checking if the last trial in this
        # was an interruption
        if not len(selected_trials):
            return np.nan, 0

        last_trial = selected_trials.iloc[-1]
        n_trials = len(selected_trials)

        if last_trial["Interrupt"] == True:
            n_trials += int(np.sum(selected_trials["RelInformativeTrialsSeen"] == last_trial["RelInformativeTrialsSeen"])) + 1

        if np.issubdtype(type(k), np.integer):
            n_bins = 1
        else:
            n_bins = len([_k for _k in k if _k <= last_trial["RelInformativeTrialsSeen"]])

        return n_trials, n_bins

    @cache
    def _p_int(self, subject: str, vocalizer: str, k: Union[int, Iterable[int]]):
        """Emperical probability of interruption

        This is the empirical probability that `subject` interrupts `vocalizer`
        after having seen `k` informative trials of that vocalizer.

        Returns
        -------
        p : float
            The probability between 0 and 1 that the subject interrupts vocalizer
            after having seen k informative trials of that vocalizer
        """
        t, n_bins = self._t_svk(subject, vocalizer, k)

        if n_bins == 0:
            return np.nan

        return (t - n_bins) / t

    @cache
    def _n_int(self, subject: str, vocalizer: str, k: Union[int, Iterable[int]]):
        """Emperical probability of interruption

        This is the empirical probability that `subject` interrupts `vocalizer`
        after having seen `k` informative trials of that vocalizer.

        Returns
        -------
        n_int, n_total : tuple of integers
            The number of interruptions of that vocalizer after seing k interupted trials
        """
        t, n_bins = self._t_svk(subject, vocalizer, k)

        if n_bins == 0:
            return 0, 0

        return t-n_bins, t

    @cache
    def p(self, subject: str, k: Union[int, Iterable[int]]):
        """Probability of interruption averaged over vocalizers, for given subject

        Implementation of (Equation 3) when dataframe in Tsvk object is limited
        to either Rewarded or Nonrewarded trials
        
        Returns
        -------
        p : float
            The probability between 0 and 1 that the subject interrupts a vocalizer
            in the selected set of trials after having seen k informative trials.
        """
        return np.nanmean([
            self._p_int(subject, v, k)
            for v in self.subject_to_vocalizers[subject]
        ])

    @cache
    def n(self, subject: str, k: Union[int, Iterable[int]]):
        """Probability of interruption averaged over vocalizers, for given subject

        Implementation of (Equation 3) when dataframe in Tsvk object is limited
        to either Rewarded or Nonrewarded trials
        
        Returns
        -------
        n : float
            The probability between 0 and 1 that the subject interrupts a vocalizer
            in the selected set of trials after having seen k informative trials.
        """
        nint = 0
        ntot = 0
        for v in self.subject_to_vocalizers[subject]:
            n1, n2 = self._n_int(subject, v, k)
            nint += n1
            ntot += n2
        
        return nint, ntot


    @cache
    def p_by_subjects(self, k: Union[int, Iterable[int]]) -> pd.DataFrame:
        """Probability of interruption averaged over vocalizers, for each subject

        This is the left hand side of of (Equation 3) for a fixed k, for each subject
        in the dataset

        Returns
        -------
        result : pd.DataFrame
            A pandas dataframe with 2 columns, "Subject" and "P_int"

            Subject : str
                A list of subject names
            P_int : float
                The avg probability of interruption for given subject after the
                kth informative trial and before the k+1th; i.e. the result
                of Equation 3 for each subject and fixed k
        """
        return pd.DataFrame({
            "Subject": self.subjects,
            "P_int": [
                self.p(subject=s, k=k)
                for s in self.subjects
            ]
        })

    def p_by_k(self) -> pd.DataFrame:
        """Estimates the probability of interruption over all subjects for each informative trial bin

        This is used to plot the probability curves in Figure 3A, 3B

        Returns
        -------
        result : pd.DataFrame
            A pandas dataframe with 3 columns, "k", "P_int" and "SE"

            k : int
                The informative trial bins k
            P_int : float
                The jackknifed estimate of the mean probability of interruption over subjects
                for each informative trial bin, from 0 to 1
            SE : float
                The standard error of the mean for the estimated interruption probability,
                computed using jackknife over subjects
        """
        means = []
        sems = []
        for k in self.k:
            p = self.p_by_subjects(k)["P_int"]
            if np.any(np.isnan(p)):
                logger.warn(f"Found a NaN element in estimate of p_by_k at k={k}")

            mean, sem = jackknife(
                p[~np.isnan(p)],
                np.mean
            )
            means.append(mean)
            sems.append(sem)
        
        return pd.DataFrame({
            "k": self.k,
            "P_int": means,
            "SE": sems
        })

    @cache
    def odds(self, subject: str, k: Union[int, Iterable[int]]) -> float:
        """Odds of interruption when p is averaged over vocalizers, for given subject

        Application of definition of odds (i.e. Odds = p / (1-p)) to interruption
        probabilities of (Equation 2)

        Handles cases where p may equal 1 or 0
        
        Returns
        -------
        odds : float
            The odds in range (-inf, inf) of interruption for subject 
        """
        p = self.p(subject, k)
        if p == 1:
            raise RuntimeError("odds() got value of p = 1, which should not be possible")
        
        if p == 0:
            # This subject never interrupted any vocalizer at k informative trials.
            # To avoid discontinuity in odds ratios, map p=0 to an estimate based on the
            # average probability of interruption of all vocalizers.
            other_probs = [
                self.p(subject, k)
                for s2 in self.subjects
                if (subject != s2) and (self.p(subject, k) != 0) and not np.isnan(self.p(subject, k))
            ]

            if not len(other_probs):
                # Return 1 / (2 * len(subjects x vocalizers))
                n_subjects_by_vocalizers = np.sum([len(v) for v in self.subject_to_vocalizers.values()])
                p = 1.0 / (2.0 * n_subjects_by_vocalizers)
                logger.debug(f"Encountered case where p=0 for all vocalizers, setting p={p}")
            else:
                p = 0.5 * np.mean(other_probs)

        return p / (1.0 - p)

    @cache
    def odds_by_subjects(self, k: Union[int, Iterable[int]]) -> pd.DataFrame:
        """Odds of interruption averaged over vocalizers, for each subject

        This is the definition of odds = p / (1 - p) applied to the
        left hand side of of (Equation 3) for a fixed k, for each subject
        in the dataset

        Returns
        -------
        result : pd.DataFrame
            A pandas dataframe with 2 columns, "Subject" and "Odds"

            Subject : str
                A list of subject names
            Odds : float
                The odds of interruption for given subject after the
                kth informative trial and before the k+1th; i.e. the result
                of Equation 3 for each subject and fixed k
        """
        return pd.DataFrame({
            "Subject": self.subjects,
            "Odds": [self.odds(subject, k) for subject in self.subjects]
        })

    @cache
    def logOR_by_subjects(self, k: Union[int, Iterable[int]]) -> pd.DataFrame:
        """Compute the logOR by subject

        This implements (Equation 4) for a given k

        Parameters
        ----------
        k : int
            The informative trial bin to compute the logOR in

        Returns
        -------
        logOR : pd.DataFrame
            A pandas dataframe with columns "Subject" and "logOR".
            logOR is the base 2 logarithm of odds ratio between interruptions
            of Nonrewarded to Rewarded vocalizers
        """
        if not self._includes_both_reward_classes:
            raise ValueError("Cannot compute logOR for a dataframe with only Rewarded or Nonrewarded trials. "
                    "This instance of Tsvk has {re_count} rewarded and {nore_count} nonrewarded trials. "
                    "Double check that you are calling `T.logOR_by_subjects` and not `T.re.logOR_by_subjects` "
                    "or `T.nore.logOR_by_subjects`".format(
                        re_count=len(self.df[self.df["StimulusClass"] == "Rewarded"]),
                        nore_count=len(self.df[self.df["StimulusClass"] == "Nonrewarded"]),
                    ))
        
        resultDict = { "Subject" : [],
                     "logORp": [],
                     "logORc": [],
                     "noreCounts": [],
                     "reCounts": []
                     }
        
        for subject in self.subjects:
            resultDict["Subject"].append(subject)
            
            noreOdds = self.nore.odds(subject, k)
            reOdds = self.re.odds(subject, k)
            resultDict["logORp"].append(np.log2(noreOdds/reOdds))
    
    
            noreNint = self.nore.n(subject, k)[0]
            if noreNint == 0:
                noreNint = 0.5
            noreNnoint = self.nore.n(subject, k)[1] - noreNint
            reNint = self.re.n(subject, k)[0]
            if reNint == 0:
                reNint = 0.5           
            reNnoint = self.re.n(subject, k)[1] - reNint
            resultDict["logORc"].append(np.log2(noreNint*reNnoint/(noreNnoint*reNint)))
            
            resultDict["noreCounts"].append(self.nore.n(subject,k))
            resultDict["reCounts"].append(self.re.n(subject,k))
            
    
        return pd.DataFrame(resultDict)
 
            
    def logOR(self, mode='average-pvalue-vocalizer'):
        """Computes the logOR over all informative trials with stats

        This is used to construct the informative trial curves in the logOR space in Figure 3A

        Returns
        -------
        logOR : pd.DataFrame
            A pandas dataframe with columns "k", "logOR", "SE", "pvalue", "tstat", and "dof"

            k : int
                The informative trial bin
            logOR : float
                The estimated log (base 2) odds ratio of interrupting Nonrewarded to interrupting
                Rewarded vocalizers, jack-knifing over subjects
            SE : float
                The standard error of the mean logOR estimated by jack-knifing over subjects
            pvalue : float
                The p-value resulting from a one-sided t-test comparing the log odds of interrupting
                NoRe to Re vocalizers over subjects
            tstat : float
                The t-statistic resulting from the one-sided t-test
            dof : int
                The degrees of freedom (n_subjects - 1) of the one-sided t-test
        """
        if mode not in ("fisher-exact", "average-logOdds-subject", "average-pvalue-vocalizer"):
            raise ValueError

        means = []
        sems = []
        pvalues = []
        tstats = []
        dofs = []
        sub_pvalues = []
        
        if mode == "fisher-exact":
            # This mode perforrms the fisher exact test for each subject based on the
            # total number of interrupts and a single fisher exact test given all the data.
            
            # This measure is biased when
            # interruption rates are low, becauase each vocalizers are under-represented
            # in the contingency matrix when their interruption rates are low (few
            # trials within each informative trial bin)
            #
            # Example:
            # v1: 4 interrupts, 1 non-interrupt, p(int|v1) = 0.8
            # v2: 0 interrupts, 1 non-interrupt, p(int|v2) = 0.0
            # v3: 0 interrupts, 1 non-interrupt, p(int|v3) = 0.0
            #
            # If we assume a uniform prior over vocalizers, then p(v_i) = 1/3
            # and p(int) = p(int|v)p(v) = (1/3) * (0.8) = 27%
            #
            # However, if we build a contingency matrix, the matrix will see this
            # probability as
            # p(int) = (4 interrupts) / (7 trials) = 57%
            #
            # This effect is most severe for low interrupting subjects, where there are
            # often 0 or 1 interruptions in a bin. So we prefer the "average" method
            # of computing logOR in an informative trial bin.

            for k in self.k:
                pvalueSubjects = np.zeros(len(self.subjects))
                tableAll = np.zeros((2,2))
                for i, subject in enumerate(self.subjects):
                    odds_ratio, ci, pvalueSubject, table = self.fisher_exact(subject, k=k, side="greater")  
                    pvalueSubjects[i] = pvalueSubject
                    tableAll = tableAll + table
            
                odds_ratio, ci_95, pvalue, se = fisher_exact(tableAll)
                mean = np.log2(odds_ratio)
                sem = (1/np.log(2))*se


                means.append(mean)
                sems.append(sem)
                pvalues.append(pvalue)
                tstats.append([])
                sub_pvalues.append(pvalueSubjects)
                dofs.append([])

        elif mode == "average-logOdds-subject":
            # This mode computes the odds and thus the log odds by on a subject basis by
            # summing all interrupted trials across all vocalizers.  
   

            for k in self.k:
                odds_ratios = np.zeros(len(self.subjects))
                for i, subject in enumerate(self.subjects):
                    odds_ratio, ci, pvalueSubject, table = self.fisher_exact(subject, k=k, side="greater")  
                    odds_ratios[i] = odds_ratio
            
                
                # This averaging the log of odds ratio
                log_odds_ratio = np.log2(odds_ratios)

                if np.any(np.isnan(log_odds_ratio)):
                    logger.warn(f"Found a NaN element in estimate of log odds ratio at k={k}")
                mean, sem = jackknife(log_odds_ratio[~np.isnan(log_odds_ratio)], np.mean)

                # One-sided t-test
                tstat, pvalue = scipy.stats.ttest_1samp(
                    log_odds_ratio[~np.isnan(log_odds_ratio)],
                    0,
                    alternative="greater",
                )
                dof = len(self.subjects) - 1

                means.append(mean)
                sems.append(sem)
                pvalues.append(pvalue)
                tstats.append(tstat)
                dofs.append(dof)
                sub_pvalues.append([])

                
        elif mode == 'average-pvalue-vocalizer':
            for k in self.k:
                log_odds_nore = np.log2(self.nore.odds_by_subjects(k)["Odds"])
                log_odds_re = np.log2(self.re.odds_by_subjects(k)["Odds"])

                # Jack-knife the expectation in Equation 4 and estimated standard error
                log_odds_ratio = log_odds_nore - log_odds_re
                if np.any(np.isnan(log_odds_ratio)):
                    logger.warn(f"Found a NaN element in estimate of log odds ratio at k={k}")
                mean, sem = jackknife(log_odds_ratio[~np.isnan(log_odds_ratio)], np.mean)

                # One-sided t-test
                tstat, pvalue = scipy.stats.ttest_rel(
                    log_odds_nore,
                    log_odds_re,
                    alternative="greater",
                )
                dof = len(self.subjects) - 1

                means.append(mean)
                sems.append(sem)
                pvalues.append(pvalue)
                tstats.append(tstat)
                dofs.append(dof)
                sub_pvalues.append([])

        return pd.DataFrame({
            "k": self.k,
            "logOR": means,
            "SE": sems,
            "pvalue": pvalues,
            "Subject Pvalues": sub_pvalues,
            "tstat": tstats,
            "dof": dofs
        })

    def fisher_exact_by_subjects(self, k: int = None, side: str = "two.sided"):
        """Perform a fisher exact test for each subject
        """
        rows = []
        for subject in self.subjects:
            odds_ratio, ci, pvalue, table = self.fisher_exact(subject, k=k, side=side)
            rows.append({
                "Subject": subject,
                "logOR": np.log2(odds_ratio),
                "95CI_low": np.log2(ci[0]),
                "95CI_high": np.log2(ci[1]),
                "pvalue": pvalue,
                "table": table
            })
        return pd.DataFrame(rows)

    @cache
    def fisher_exact(self, subject: str, k: int = None, side: str = "two.sided"):
        """Perform a the fisher exact test on rewarded and non-rewarded vocalizers

        If k is specified, perform it for only that informative trial bin, which
        includes all trials with k informative trials seen.

        Otherwise, perform it on the entire dataset.

        Returns
        -------
        odds_ratio : float
            Estimate of the odds ratio
        ci_95 : tuple[float, float]
            95% confidence interval bounds
        p_value : float
            P value for significance testing
        """
        if side not in ("greater", "less", "two.sided"):
            raise ValueError("side parameter must be one of 'greater', 'less', or 'two.sided'")
        # If k is None
        if k is None:
            re_df = self.re.df
            nore_df = self.nore.df
        else:
            if np.issubdtype(type(k), np.integer):
                re_df = self.re.df[self.re.df["RelInformativeTrialsSeen"] == k]
                nore_df = self.nore.df[self.nore.df["RelInformativeTrialsSeen"] == k]
            else:
                re_df = self.re.df[self.re.df["RelInformativeTrialsSeen"].isin(k)]
                nore_df = self.nore.df[self.nore.df["RelInformativeTrialsSeen"].isin(k)]

        re_df = re_df[re_df["Subject"] == subject]
        nore_df = nore_df[nore_df["Subject"] == subject]

        table = np.array([
            [
                np.sum(nore_df["Interrupt"] == True),
                np.sum(nore_df["Interrupt"] == False),
            ],
            [
                np.sum(re_df["Interrupt"] == True),
                np.sum(re_df["Interrupt"] == False),
            ]
        ])

        odds_ratio, ci_95, p_value, se = fisher_exact(table)
        return odds_ratio, ci_95, p_value, table
