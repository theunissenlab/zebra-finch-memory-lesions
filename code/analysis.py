import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity

from load_data import load_ephys_stimulus
from stats import jackknife
from utils import generate_synthetic_poisson_spike_times


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
    bandwidth (float): bandwidth of gaussian kernel
    return_rate (bool): return kde as estimate of firing rate if True. If False, returns
        kde as a probability distribution
    estimate_std_error (bool, default False): use jackknife procedure
        to estimate standard error of the kde. Is much slower, so only use if you need it.

    Returns
    =======
    t_arr: numpy array of time in seconds spanning time_range
    spike_density: the psth in normalized units (if return_rate was False) or in Hz
    """
    n_samples = int(np.round((time_range[1] - time_range[0]) * psth_sampling_rate))
    t_kde = (np.arange(n_samples) / psth_sampling_rate) + time_range[0]

    if estimate_std_error:
        log_density, log_density_se = jackknife(spike_times, _log_density, t_kde=t_kde, bandwidth=bandwidth)
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
    def _max_kde(spike_times, stim_duration=None, return_rate=True):
        t_arr, result, (ci_low, ci_high) = fit_kde(spike_times, (0, stim_duration), return_rate=True)
        return np.max(result)

    @staticmethod
    def max_response(spike_times, stim_duration, compute_std_error=False):
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


__all__ = [
    "ResponseStrength",
    "fit_kde",
    "get_stimulus_timeseries",
]
