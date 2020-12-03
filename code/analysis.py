import logging

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn.neighbors import KernelDensity

from load_data import load_ephys_stimulus
from stats import false_discovery, jackknife
from utils import clean_spike_times, generate_synthetic_poisson_spike_times


logger = logging.getLogger(__name__)


def get_stimulus_timeseries(stimulus_file, time_range):
    """Create time and amplitude arrays for a stimulus file over a given time range

    The stimulus onset will be aligned to time=0.0

    Params
    ======
    stimulus_file (string or pandas.core.series.Series): name of stimulus file to load OR
        row of already loaded stimulus dataframe
    time_range (tuple): tuple of start time and stop time of data to return in seconds relative
        to trial onset
    """
    if isinstance(stimulus_file, pd.core.series.Series):
        fs = stimulus_file["sampling_rate"]
        stim = stimulus_file["stim"]
        duration = stimulus_file["duration"]
    else:
        wav = load_ephys_stimulus(stimulus_file)
        fs = wav.rate
        stim = wav.data[:, 0]
        duration = len(stim) / fs

    if time_range[0] > 0.0:
        raise ValueError("I was too lazy to implement this case")

    full_duration = time_range[1] - time_range[0]

    pad_before = np.zeros(int(fs * -time_range[0]))
    n_pad_after = int(fs * full_duration) - len(pad_before) - stim.size
    if n_pad_after > 0:
        pad_after = np.zeros(n_pad_after)
        stim = np.concatenate((pad_before, stim, pad_after))
    else:
        stim = stim[:n_pad_after - 1]
        stim = np.concatenate((pad_before, stim))
    t_arr = np.linspace(time_range[0], full_duration + time_range[0] - 1 / fs, stim.size)

    return t_arr, stim


def _log_density(spike_times, t_kde=None, bandwidth=0.03):
    """Use sklearn to estimate log density of spike_times

    Params
    ======
    t_kde (np.ndarray): array of time points to sample log density at
    bandwidth (float): bandwidth of gaussian kernel for kde estimate
    """
    spikes = np.zeros((0, 1))
    if not len(spike_times):
        spikes = np.zeros((0, 1))
    else:
        spikes = np.concatenate(spike_times)
        spikes = spikes.reshape(-1, 1)

    if not len(spikes):
        return np.full(t_kde.shape, -np.inf)

    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth, rtol=1e-2).fit(spikes)
    return kde.score_samples(t_kde[:, None])


def fit_kde(
        spike_times,
        time_range,
        psth_sampling_rate=1000,
        bandwidth=0.03,
        return_rate=True,
        estimate_std_error=False
    ):
    """Estimates a psth with a kernel density estimate using a gaussian kernel

    Params
    ======
    spike_times (list of N lists): spike times (in seconds from trial onset) for N trials
    time_range (tuple): tuple of start time and stop time of psth in seconds relative to
        trial onset to generate psth from
    psth_sampling_rate (float): sampling rate of generated psth
    bandwidth (float, default 0.03): bandwidth of gaussian kernel.
        If "auto", uses 0.1 divided by number of trials. Otherwise, uses float value provided
    return_rate (bool): return kde as estimate of firing rate if True. If False, returns
        kde as a probability distribution
    estimate_std_error (bool, default False): use jackknife procedure
        to estimate standard error of the kde. Is much slower, so only use if you need it.

    Returns
    =======
    t_arr: numpy array of time in seconds spanning time_range
    spike_density: the psth in normalized units (if return_rate was False) or in Hz
    """
    if bandwidth == "auto":
        bandwidth = 0.1 / len(spike_times)

    n_samples = int(np.round((time_range[1] - time_range[0]) * psth_sampling_rate))
    t_kde = (np.arange(n_samples) / psth_sampling_rate) + time_range[0]

    if estimate_std_error:
        logging.warn("The jackknife method for the kde fit doesn't seem to work very well...")
        log_density, log_density_se = jackknife(spike_times, _log_density, t_kde=t_kde, bandwidth=bandwidth)
        spike_probability_density = np.exp(log_density)
        ci_95_lower = np.exp(
            log_density - 2 * log_density_se
        )

        ci_95_upper = np.exp(
            log_density + 2 * log_density_se
        )
    else:
        log_density = _log_density(spike_times, t_kde=t_kde, bandwidth=bandwidth)
        spike_probability_density = np.exp(log_density)
        ci_95_lower = None
        ci_95_upper = None

    spikes = np.concatenate(spike_times)

    if return_rate:
        n_spikes = len(spikes)
        n_trials = len(spike_times)
        spike_density = spike_probability_density * n_spikes / n_trials
        if ci_95_lower is not None:
            ci_95_lower *= n_spikes / n_trials
            ci_95_upper *= n_spikes / n_trials
    else:
        spoike_density = spike_probability_density

    return t_kde, spike_density, (ci_95_lower, ci_95_upper)


class ResponseStrength:
    """Organization object for functions that compute response strengths

    TODO (kevin): should baselines be calculated here as well?
    """

    @staticmethod
    def count_spikes_in_time_window(spike_times, time_window=(0, 0.5)):
        filtered_spike_rates = [
            len(s[(s >= time_window[0]) & (s < time_window[1])]) / (time_window[1] - time_window[0])
            for s in spike_times
        ]
        return filtered_spike_rates

    @staticmethod
    def rate(spike_times, time_window=(0, 0.5)):
        """Calculates a response strength based on a time window

        Calculates this average across all trials (i.e with no grouping by stim)

        Params
        ======
        spike_times: list of lists of spike times
        time_window (tuple): start and end times of time window to evaluate firing rate in

        Returns
        =======
        Returns a mean and standard deviation of the estimated firing rate during the specified window
        """
        if not len(spike_times):
            return None, None

        filtered_spike_rates = ResponseStrength.count_spikes_in_time_window(
            spike_times, time_window=time_window
        )
        return jackknife(filtered_spike_rates, np.mean, parallel=False)

    @staticmethod
    def _max_kde(
            spike_times,
            stim_duration=None,
            return_rate=True,
            normalize_to_trials=None,
            ):
        t_arr, result, (ci_low, ci_high) = fit_kde(spike_times, (0, stim_duration), return_rate=True)
        return np.max(result)

    @staticmethod
    def max_response(spike_times, stim_duration, normalized_trial_count=None, compute_std_error=False):
        """Calculates a response strength based on max kde

        Note that for low firing rates, the max kde is expected to be
        much less for short stims than long stims. So you should compute some baseline values
        to subtract off of this.

        Params
        ======
        spike_times: list of lists of spike times
        stim_duration: window in which to determine the max kde (should be no larger than the shortest stim...)

        Returns
        =======
        Returns a max kde and jackknifed standard error
        """
        if normalized_trial_count is None:
            pass
        else:
            _random_indexes = np.random.choice(
                np.arange(len(spike_times)),
                normalized_trial_count
            )
            spike_times = np.array(spike_times)[_random_indexes]

        if compute_std_error:
            mean, se = jackknife(spike_times, ResponseStrength._max_kde, stim_duration=stim_duration, return_rate=True)
        else:
            mean = ResponseStrength._max_kde(spike_times, stim_duration, return_rate=True)
            se= None

        return mean, se

    @staticmethod
    def max_response_baseline(spike_times, stim_duration, n_trials, iterations=40):
        """Estimate the expected max value of n_trials of a poisson spiking process where
        the mean firing rate is estimated from the time_window of the given spike_times

        Returns a float which is the expected value of the max kde of n_trials

        Maybe Deprecate this if I'm just gonna manually choose max_response duration parameter
        to be the same for all calls in the test
        """
        baseline_rate, _ = ResponseStrength.rate(spike_times, time_window=(0, stim_duration))

        _maxes = []
        for _ in range(iterations):
            spikes = generate_synthetic_poisson_spike_times(baseline_rate, stim_duration, n_trials)
            t, kde, _ = fit_kde(spikes, (0, stim_duration), return_rate=True)
            response, _ = ResponseStrength.max_response(spikes, stim_duration=stim_duration)
            _maxes.append(response)

        return np.mean(_maxes)


class Auditory:

    @staticmethod
    def test_auditory_by_mean_rate(unit_df, delta_t=0.5):
        """Perform a t-test to determine if the onset response is significant
        """
        spike_times = clean_spike_times(unit_df["spike_times"])
        if not len(spike_times):
            return 1.0

        baseline_spike_counts = ResponseStrength.count_spikes_in_time_window(spike_times, time_window=(-delta_t, 0.0))
        stim_spike_counts = ResponseStrength.count_spikes_in_time_window(spike_times, time_window=(0.0, delta_t))

        return ttest_ind(baseline_spike_counts, stim_spike_counts).pvalue

    @staticmethod
    def is_auditory_by_mean_rate(unit_df, alpha=0.05, delta_t=0.5):
        """Return true if the spiking data indicates a significant onset response

        Averaged across all stimuli
        """
        return Auditory.test_auditory_by_mean_rate(unit_df, delta_t=delta_t) <= alpha

    @staticmethod
    def _test_stim_auditory_by_any_rate(stim_row, baseline_counts, alpha=None, delta_t=0.5):
        """Test whether spike counts for given stimulus differs from baseline in any time window

        Breaks up the stimulus into windows of time delta_t, and performs a t-test for spike
        counts in all windows against a baseline set of spike counts. Either returns the minimum
        p-value (multiplied by number of windows to adjust for multiple comparisons), or
        the maximum significant p-value after a Benjamini-Hochberg correction for multiple
        comparisons if a significant threshold is predefined.

        Example
        -------
         Window 1  Window 2  Window 3  Window 4
        [ p=0.10 ][ p=0.01 ][ p=0.02 ][ p=0.05 ]

        * If alpha is None, will return the smallest p-value multipled by number of windows (i.e.
          p = 0.01 * 4 = 0.04).
        * If alpha is passed in as 0.05, Windows 2 and 3 will be found to be significant after
          multiple comparisons correction, so the max significant p-value will be returned (i.e.
          p = 0.02)

        Params
        ------
        stim_row : pandas.core.series.Series
            The dataframe row for the select stim. Should have an element "stim_duration",
            a float indicating the stim duration in seconds, and "spike_times", a list of
            lists of spike arrival times in seconds relative to trial onset.
        baseline_counts : list
            A list of counts of spike arrival times to compare spike counts for the given stim to
        alpha : float, default=None
            Significance threshold for multiple comparisons correction. If not specified, will
            return the smallest p-value across all windows multiplied by the number of time windows
        delta_t : float, default=0.5
            Size of time windows to test

        Returns
        -------
        A p-value describing the outcome of the test whether spike counts in any time window
        for this stimulus differ from baseline
        """
        # Test all delta_t time windows of stim
        windows = int((stim_row["stim_duration"] // delta_t) + 1)

        stim_pvalues = np.ones(windows)
        for i, t_window in enumerate(np.arange(windows) * delta_t):
            stim_spike_counts = ResponseStrength.count_spikes_in_time_window(
                stim_row["spike_times"],
                time_window=(t_window, t_window + delta_t)
            )
            stim_pvalues[i] = ttest_ind(baseline_counts, stim_spike_counts).pvalue

        # Multiple comparisons correction by multiplying the best pvalue by the number of
        # comparisons (number of time windows)?
        # I find that this is often too strict. Maybe I should do the max significant p-value
        # under false discovery correction, with no multiple
        if alpha is not None:
            significance_mask = false_discovery(stim_pvalues, alpha=alpha)
            if np.any(significance_mask):
                return(np.max(stim_pvalues[significance_mask]))

        return np.min(stim_pvalues) * windows

    @staticmethod
    def _count_spikes_in_all_preceding_time_windows(unit_df, delta_t=0.5):
        """Breaks up the silent period before each stim into windows of fixed size and reports spike counts
        """

    @staticmethod
    def test_auditory_by_any_rate(unit_df, alpha=None, delta_t=0.5):
        spike_times = clean_spike_times(unit_df["spike_times"])
        if not len(spike_times):
            return []

        # This code uses the -delta_t time period as the baseline.
        # The commented code below was an attempt to cut the silent period
        # Preceding the stims into windows of delta_t (i.e. get more baseline
        # windows)
        baseline_counts = ResponseStrength.count_spikes_in_time_window(
            spike_times,
            time_window=(-delta_t, 0)
        )

        # Get spike times in all baseline time windows for every row
        # Each row may have a different amount of
        # windows = -np.arange(1, 1 + np.floor(1.0 / delta_t)) * delta_t
#         if "preceding_silence_duration" not in unit_df.columns:
#             windows = -np.arange(1, 1 + np.floor(1.0 / delta_t)) * delta_t
#             baseline_counts = np.concatenate([
#                 ResponseStrength.count_spikes_in_time_window(
#                     spike_times,
#                     time_window=(window, window + delta_t)
#                 ) for window in windows
#             ])
#         else:
#             baseline_counts = []
#             for i in range(len(unit_df)):
#                 row = unit_df.iloc[i]
#                 windows = -np.arange(1, 1 + np.floor(row["preceding_silence_duration"] / delta_t)) * delta_t
#                 new_counts = np.concatenate([
#                     ResponseStrength.count_spikes_in_time_window(
#                         [row["spike_times"]],
#                         time_window=(window, window + delta_t)
#                     ) for window in windows
#                 ])
#                 baseline_counts = np.concatenate([baseline_counts, new_counts])

        pvalues = []
        for i in range(len(unit_df)):
            row = unit_df.iloc[i]
            pvalues.append(Auditory._test_stim_auditory_by_any_rate(
                row,
                baseline_counts,
                alpha=alpha,
                delta_t=delta_t
            ))

        return pvalues

    @staticmethod
    def is_auditory_by_any_rate(unit_df, alpha=0.05, delta_t=0.5):
        """Return true if the spiking data indicates a significant auditory response to any stim
        """
        pvalues = Auditory.test_auditory_by_any_rate(unit_df, delta_t=delta_t, alpha=alpha)
        return np.any(false_discovery(pvalues, alpha=alpha))

    @staticmethod
    def fraction_auditory_by_any_rate(unit_df, alpha=0.05, delta_t=0.5):
        """Return fraction (0.0 to 1.0) of stims with a significant auditory at any time
        """
        pvalues = Auditory.test_auditory_by_any_rate(unit_df, delta_t=delta_t)
        return np.mean(false_discovery(pvalues, alpha=alpha))

    @staticmethod
    def selectivity(unit_df, mode="rate", **kwargs):
        """Calculate selectivity for all stims

        The selectivity to each stim is calculated as the ratio
        of the response strength to that stim divided by the mean response
        strength across all stims
        """
        if mode not in ("rate", "max_response", "n_auditory"):
            raise ValueError("mode must be either 'rate', 'max_response', or 'n_auditory'")

        spike_times = [s for s in list(unit_df["spike_times"]) if len(s)]
        if not len(spike_times):
            return np.array([])

        if mode == "rate":
            if "time_window" not in kwargs:
                kwargs["time_window"] = (0, 0.5)
            spike_times = np.concatenate(spike_times)
            mean_response = ResponseStrength.rate(spike_times, time_window=(0.0, 0.5))[0]

            stim_responses = np.array([
                ResponseStrength.rate(
                    unit_df.iloc[i]["spike_times"], **kwargs)[0]
                    if len(unit_df.iloc[i]["spike_times"])
                    else 0
                for i in range(len(unit_df))
            ])
        elif mode is "max_response":
            if "stim_duration" not in kwargs:
                kwargs["stim_duration"] = np.min(unit_df["stim_duration"])

            stim_responses = np.array([
                # We used to have a check that the length of spike_times was not zero...
                ResponseStrength.max_response(unit_df.iloc[i]["spike_times"], **kwargs)[0]
                for i in range(len(unit_df))
            ])
            mean_response = np.mean(stim_responses)

        return stim_responses / mean_response

    @staticmethod
    def selectivity_index(unit_df, mode="rate", **kwargs):
        """Return selectivity index for the given stimuli

        The selectivity index (SI) is max response to a single stim divided
        by the mean response over all stims.

        Params
        ======
        unit_df : pd.DataFrame
            rows of stimulus response data
        mode : string (default="rate")
            either "rate", "max_response", or "n_auditory"
            kwargs specified are passed into the relevant function
        kwargs : dict
            arguments passed into the response functions to test for auditoryness
            for rate, include a "time_window" tuple
            for max_response, include a "duration"
            for n_auditory, include an "alpha" value for significance tests and "delta_t"
                for time window widths
        """
        if mode not in ("rate", "max_response", "n_auditory"):
            raise ValueError("mode must be either 'rate', 'max_response', or 'n_auditory'")

        if mode == "n_auditory":
            return Auditory.fraction_auditory_by_any_rate(unit_df, **kwargs)
        else:
            selectivity = Auditory.selectivity(unit_df, mode=mode, **kwargs)
            return np.max(selectivity)



__all__ = [
    "ResponseStrength",
    "fit_kde",
    "get_stimulus_timeseries",
]
