"""Create a datf not os.path.exists(target_file):
    unit_summary.to_pickle(target_file)aframe with auditory information for each unit
"""

import os
import sys
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(PROJECT_ROOT, "code"))

import numpy as np
import pandas as pd
import tqdm
from pandarallel import pandarallel

import groups
from analysis import (
    Auditory,
    ResponseStrength as RS
)
from load_data import (
    load_spike_data,
    load_stim_data,
    load_rendition_data,
)
from utils import clean_spike_times, yes_no


pandarallel.initialize(progress_bar=False)

import warnings
warnings.simplefilter('error', RuntimeWarning)


if __name__ == "__main__":
    output_path = os.path.join(PROJECT_ROOT, "data", "ephys", "AuditoryUnitTable.pkl")
    if os.path.exists(output_path):
        if not yes_no("Dataframe at {} exists.\nOverwrite? [y/n]:".format(output_path)):
            print("Aborting.")
            sys.exit(0)

    print("Loading data...")
    ephys_data = load_spike_data()
    stim_data = load_stim_data()
    rendition_data = load_rendition_data()
    print("Data loaded.")

    # Compute auditory flags for different methods
    columns = {}
    print("Calculating firing rate of all units")
    columns["fr"] = ephys_data.groupby("unit_id").apply(
        lambda unit_df: RS.rate(clean_spike_times(unit_df["spike_times"]))[0]
    )
    print("Finding units with auditory responses in any 1s window")
    columns["auditory_in_any_1s"] = ephys_data.groupby("unit_id").parallel_apply(
        lambda unit_df: Auditory.is_auditory_by_any_rate(unit_df, alpha=0.01, delta_t=1.0),
    )
    print("Finding units with auditory responses in any 200ms window")
    columns["auditory_in_any_200ms"] = ephys_data.groupby("unit_id").parallel_apply(
        lambda unit_df: Auditory.is_auditory_by_any_rate(unit_df, alpha=0.01, delta_t=0.2)
    )
    print("Finding units with audiotory responses in first 500ms window")
    columns["auditory_onset_500ms"] = ephys_data.groupby("unit_id").parallel_apply(
        lambda unit_df: Auditory.is_auditory_by_mean_rate(unit_df, alpha=0.05),
    )

    # Create a table for the rate data
    rate_auditory_table = pd.DataFrame(
        columns,
        index=ephys_data.groupby("unit_id").grouper.levels[0]
    ).reset_index()
    auditory_ephys_data = ephys_data.merge(rate_auditory_table, on="unit_id")

    # Collapse data into one row per unit with additional unit information
    unit_summary = auditory_ephys_data.groupby("unit_id").agg({
        "subject": lambda x: x[0],
        "electrode": lambda x: x[0],
        "unit": lambda x: x[0],
        "site": lambda x: x[0],
        "snr": lambda x: x[0],
        "n_trials": "mean",
        "isi_violation_pct": lambda x: x[0],
        "fr": lambda x: x[0],
        "auditory_onset_500ms": lambda x: x[0],
        "auditory_in_any_200ms": lambda x: x[0],
        "auditory_in_any_1s": lambda x: x[0],
    })

    # Compute firing rates just to all stimuli and just to songs, dcs, ripples
    print("Calculating firing rate for songs, dcs, and ripples")
    splits = [groups.all, groups.songs, groups.dcs, groups.ripples]
    rate_data = []
    for i in tqdm.tqdm(range(len(unit_summary)), total=len(unit_summary)):
        row = unit_summary.iloc[i]
        _rendition_unit_data = rendition_data.query("unit_id == '{}'".format(row.name))
        result = {}
        for split in splits:
            fr,  _ = RS.rate(
                clean_spike_times(split(_rendition_unit_data)["spike_times"]),
                time_window=(0, 0.5)
            )
            result["{}_fr".format(split.__name__)] = fr
        rate_data.append(result)

    for split in splits:
        key = "{}_fr".format(split.__name__)
        unit_summary[key] = pd.Series([x[key] for x in rate_data], index=unit_summary.index)

    unit_summary.to_pickle(output_path)
