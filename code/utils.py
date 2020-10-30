import numpy as np


def next_poisson_dt(rate):
    """Generate a spike time interval based on firing rate"""
    return -np.log(np.random.random()) / rate


def generate_synthetic_poisson_spike_times(firing_rate, duration, n_trials):
    trials = []
    for i in range(n_trials):
        trial = []
        curr_t = 0
        while curr_t < duration:
            dt = next_poisson_dt(firing_rate)
            if curr_t + dt > duration:
                trials.append(trial)
                break
            else:
                trial.append(curr_t + dt)
                curr_t += dt
    return np.array([np.array(trial) for trial in trials])


def clean_spike_times(spike_times):
    """Collapse a list of lists of spike times while avoiding errors from empty datasets
    """
    if not len(spike_times):
        return np.array([])
    else:
        arrs = [sts for sts in spike_times if len(sts)]
        if not len(arrs):
            return np.array([])
        else:
            return np.concatenate(arrs)


class Colors:
    base_rewarded = "#0AA5D8"
    base_nonrewarded = "#C62533"
    song = "#F75F23"
    shuffled_song = "#FBA989"
    dc = "#3CAEA3"
    ripple = "#AAAB57"



__all__ = [
    "generate_synthetic_poisson_spike_times",
]
