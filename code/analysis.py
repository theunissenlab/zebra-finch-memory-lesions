import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity

from load_data import load_ephys_stimulus


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


def fit_kde(spike_times, time_range, psth_sampling_rate=1000, bandwidth=0.03, return_rate=True):
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

    Returns
    =======
    t_arr: numpy array of time in seconds spanning time_range
    spike_density: the psth in normalized units (if return_rate was False) or in Hz
    """
    if not len(spike_times):
        spikes = np.zeros((0, 1))
    else:
        spikes = np.concatenate(spike_times).reshape(-1, 1)

    n_samples = int(np.round((time_range[1] - time_range[0]) * psth_sampling_rate))
    t_kde = (np.arange(n_samples) / psth_sampling_rate) + time_range[0]
    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth, rtol=1e-2).fit(spikes)
    spike_density = np.exp(kde.score_samples(t_kde[:, None]))
    spike_density = spike_density / np.sum(spike_density)
    if return_rate:
        n_spikes = len(spikes)
        n_trials = len(spike_times)
        spike_density = spike_density * n_spikes * psth_sampling_rate / n_trials

    return t_kde, spike_density
