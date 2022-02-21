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

from .stats import jackknife


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
        """
        if isinstance(k, np.integer):
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

        if not len(selected_trials):
            # Estimate probabilty of interruption
            # let X = number of interruptions since the last informative trial
            # Return 2X (so that p = (X - 0.5) / X
            selected_trials = self.df[
                (self.df["Subject"] == subject) &
                (self.df["StimulusVocalizerId"] == vocalizer)
            ]
            last_k = np.max(selected_trials["RelInformativeTrialsSeen"])
            X = np.sum(selected_trials["RelInformativeTrialsSeen"] == last_k)
            logger.debug(f"Encountered no trials for subject={subject} | vocalizer={vocalizer} | k={k}; setting t_svk to {2*X}")
            return np.nan
            return 2 * X

        return len(selected_trials)

    @cache
    def _p_int(self, subject: str, vocalizer: str, k: Union[int, Iterable[int]]):
        """Emperical probability of interruption

        Implementation of (Equation 2)

        This is the empirical probability that `subject` interrupts `vocalizer`
        after having seen `k` informative trials of that vocalizer.

        Returns
        -------
        p : float
            The probability between 0 and 1 that the subject interrupts vocalizer
            after having seen k informative trials of that vocalizer
        """
        t = self._t_svk(subject, vocalizer, k)
        return (t - 1) / t

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
        """Odds of interruption averaged over vocalizers, for given subject

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

        return pd.DataFrame({
            "Subject": self.subjects,
            "logOR": [
                np.log2(self.nore.odds(subject, k)) - np.log2(self.re.odds(subject, k))
                for subject in self.subjects
            ]
        })
            
    def logOR(self):
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
        means = []
        sems = []
        pvalues = []
        tstats = []

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
                log_odds_re
            )
            dof = len(self.subjects) - 1

            # One-sided t-test p value conversion
            pvalue = pvalue / 2 if tstat > 0 else (1 - (pvalue / 2))

            means.append(mean)
            sems.append(sem)
            pvalues.append(pvalue)
            tstats.append(tstat)

        return pd.DataFrame({
            "k": self.k,
            "logOR": means,
            "SE": sems,
            "pvalue": pvalues,
            "tstat": tstats,
            "dof": dof
        })

