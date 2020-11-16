import os
import multiprocessing
import sys
if not os.path.basename(os.getcwd()) == "zebra-finch-memory-lesions":
    raise Exception("Run this script from the top level of the project")

sys.path.append(os.path.join(os.getcwd(), "code"))

from soundsig.sound import spectrogram, plot_spectrogram

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import groups
from analysis import (
    fit_kde,
    get_stimulus_timeseries
)
from coherence import bootstrap_coherence
from load_data import (
    load_auditory_info,
    load_spike_data,
    load_stim_data,
    load_ephys_stimulus,
    load_rendition_data,
    load_unit_waveforms
)
from plotting import plot_raster, fig_grid
from utils import Colors, clean_spike_times, generate_synthetic_poisson_spike_times


def plot_coherence_result(coherence_result, ax=None, color="Black"):
    if ax is None:
        ax = plt.gca()
    ax.plot(
        np.fft.fftshift(coherence_result["freqs"]),
        np.fft.fftshift(coherence_result["coherence"]),
        color=color,
        linestyle="-"
    )
    plt.fill_between(
        np.fft.fftshift(coherence_result["freqs"]),
        np.fft.fftshift(coherence_result["coherence_bounds"][0]),
        np.fft.fftshift(coherence_result["coherence_bounds"][1]),
        color=color,
        alpha=0.3,
    )
    
def _coherence_integrator(coherence_dict, significant_only=True):
    delta_f = coherence_dict["freqs"][1] - coherence_dict["freqs"][0]
    
    selector = coherence_dict["freqs"] >= 0
    if significant_only:
        selector = selector & (coherence_dict["coherence_bounds"][0] > 0)
    freqs = coherence_dict["freqs"][selector]
    return freqs, np.cumsum(-np.log2(1 - coherence_dict["coherence"])[selector]) * delta_f


def information(coherence_dict, min_freq=None, max_freq=None, significant_only=True):
    delta_f = coherence_dict["freqs"][1] - coherence_dict["freqs"][0]

    selector = coherence_dict["freqs"] >= 0
    if min_freq:
        selector = selector & (coherence_dict["freqs"] >= min_freq)
    if max_freq:
        selector = selector & (coherence_dict["freqs"] <= max_freq)
    if significant_only:
        selector = selector & (coherence_dict["coherence_bounds"][0] > 0)

    info = np.sum(-np.log2(1 - coherence_dict["coherence"])[selector]) * delta_f
    return info


def plot_coherence_and_information(coherence_result):
    plot_coherence_result(coherence_result)
    plt.ylim(0, 1)
    plt.xlim(0, 80)
    plt.twinx()
    plt.plot(*_coherence_integrator(coherence_result, significant_only=True))
    plt.ylim(0, 2 * plt.ylim()[1])
    plt.ylabel("Information (bits)")
    plt.text(0.9, 0.9, "I={:.2f} Bits".format(information(coherence_result, significant_only=True)),
            transform=plt.gca().transAxes,
            horizontalalignment="right",
            verticalalignment="top")
    return coherence_result
    

def compare_coherence_and_information(df, *splits, min_trials=6, figsize=(4, 2)):
    fig, axes = fig_grid(len(splits), ax_size=figsize, max_columns=5)
    
    unit_id = df.iloc[0]["unit_id"]

    results = {}
    for i, split in enumerate(splits):
        ax = axes[i]
        ax.set_title(split.__name__)
        
        split_df = split(df)
        coherence_result = bootstrap_coherence(
            split_df,
            buffer=0.2,
            min_trials=min_trials,
            iters=100,
            parallel=False
        )
        print("Computed {}/{} for {}".format(i + 1, len(splits), unit_id))
        coherence_result = plot_coherence_and_information(coherence_result)
        results[split.__name__] = coherence_result["coherence"]
        results["{}_lower".format(split.__name__)] = coherence_result["coherence_bounds"][0]
        results["{}_upper".format(split.__name__)] = coherence_result["coherence_bounds"][1]
        results["freqs"] = coherence_result["freqs"]
    results["unit_id"] = unit_id

    return fig, results


def run(unit_df):
    unit_id = unit_df.iloc[0]["unit_id"]
    splits = [
        groups.familiar_songs, 
        groups.unfamiliar_songs, 
        groups.rewarded_songs,
        groups.nonrewarded_songs,
        groups.familiar_test_songs,
        groups.familiar_nontest_songs,
        groups.familiar_dcs,
        groups.unfamiliar_dcs,
        groups.rewarded_dcs,
        groups.nonrewarded_dcs,
        groups.familiar_test_dcs,
        groups.familiar_nontest_dcs,
        groups.ripples,
        groups.songs,
        groups.dcs,
    ]

    print("RUNNING {}".format(unit_id))
    
    fig, results = compare_coherence_and_information(unit_df, *splits)
    fig.savefig("/auto/fhome/kevin/Projects/zebra-finch-memory-lesions/data/ephys/coherence_plots/{}.svg".format(unit_id), format="svg")
    plt.close(fig)
    
    return fig, results


if __name__ == "__main__":
    ephys_data = load_spike_data()
    stim_data = load_stim_data()
    auditory_info = load_auditory_info()
    rendition_data = load_rendition_data()
    print("DATA LOADED")

    auditory_unit_ids = auditory_info[
         (auditory_info["auditory_by_rate"] | auditory_info["auditory_by_any"])   
    ]["unit_id"]

    unit_summary = ephys_data.merge(auditory_info, on="unit_id").groupby("unit_id").agg({
        "subject": lambda x: x[0],
        "electrode": lambda x: x[0],
        "unit": lambda x: x[0],
        "site": lambda x: x[0],
        "snr": lambda x: x[0],
        "spike_times": lambda x: list(clean_spike_times(x)),
        "n_trials": "mean",
        "isi_violation_pct": lambda x: x[0],
        "fr": lambda x: x[0],
        "auditory_by_rate": lambda x: x[0],
        "auditory_by_any": lambda x: x[0]})

    units = unit_summary[
        (unit_summary["snr"] > 5) &
        (unit_summary["isi_violation_pct"] < 0.001) &
        (unit_summary["n_trials"] > 8) &
        ((unit_summary["auditory_by_rate"] == True) | (unit_summary["auditory_by_any"] == True))
    ]

    print("UNITS SUMMARIZED")

    # Run coherneces on rendition dataframe
    unit_dfs = [rendition_data.query(
        "unit_id == '{}'".format(unit_id)
    ) for unit_id in units.index]

    print("RUNNING COHERENCE ON AUDITORY UNIT RENDITIONS")
    with multiprocessing.Pool(7) as p:
        results = p.map(run, unit_dfs)
    pd.DataFrame([x[1] for x in results]).to_pickle("/auto/fhome/kevin/Projects/zebra-finch-memory-lesions/data/ephys/RenditionCoherences.pkl")


