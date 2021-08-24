"""Expand the data from each stimulus playback into seperate responses to each rendition

Each song playback included 3 renditions each, DCs had 3, and ripples had 1.

Save the output into a new data file
"""

import os
import sys
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(PROJECT_ROOT, "code"))

import tqdm

import numpy as np
import pandas as pd

from load_data import (
    load_spike_data,
    load_stim_data,
    load_ephys_stimulus,
)
from utils import yes_no

DATA_DIR = os.path.join(PROJECT_ROOT, "data")


def create_split_dataframe(df, stim_df, stimulus_file_to_rendition_bounds):
    output_rows = []
    for row_idx in tqdm.tqdm(range(len(df))):
        row = df.iloc[row_idx]
        stimulus_file = row["stimulus_file"]

        stim_bounds = stimulus_file_to_rendition_bounds[stimulus_file]

        prev_rendition_stop_time = None
        for rendition_idx, (rendition_start_time, rendition_stop_time) in enumerate(stim_bounds):
            new_row = dict(row)
            new_row["spike_times"] = row["spike_times"] - rendition_start_time
            new_row["stim_duration"] = rendition_stop_time - rendition_start_time
            new_row["rendition_idx"] = rendition_idx

            # Since there is no guarantee of a silence period duration before each rendition
            # We check here if there is enough time preceding this rendition to consider silence;
            # We also provide a timestamp of the trial onset.

            # A 1s buffer was defined when the ephys trial data was first exported
            # Use this for the first renditions always
            if prev_rendition_stop_time is None:
                new_row["preceding_silence_duration"] = 1.0
            else:
                new_row["preceding_silence_duration"] = (
                    rendition_start_time
                    - prev_rendition_stop_time
                )

            new_row["trial_onset"] = -rendition_start_time
            prev_rendition_stop_time = rendition_stop_time

            output_rows.append(new_row)

    return pd.DataFrame(output_rows)


if __name__ == "__main__":
    output_path = os.path.join(DATA_DIR, "ephys", "UnitRenditionData.pkl")

    if os.path.exists(output_path):
        if not yes_no("Renditon dataframe {} exists.\nOverwrite? [y/n]:".format(output_path)):
            print("Aborting.")
            sys.exit(0)

    print("Loading pickle files")
    ephys_data = load_spike_data()
    stim_data = load_stim_data()

    # Cache the onset and offset times for each rendition by stimuluse_file
    _stim_rendition_bounds = {}
    for stim_idx in range(len(stim_data)):
        stim_row = stim_data.iloc[stim_idx]
        stim_bounds = [
            stim_row["t_stim"][list(idx_bounds)] for idx_bounds in stim_row["rendition_indexes"]
        ]
        _stim_rendition_bounds[stim_row["stimulus_file"]] = stim_bounds

    import time
    run_start = time.time()
    ephys_rendition_data = create_split_dataframe(
        ephys_data,
        stim_data,
        _stim_rendition_bounds
    )

    ephys_rendition_data.to_pickle(output_path)

    print("Created dataframe split by rendition "
          "in {:.2f} seconds\n{}".format(time.time() - run_start, output_path))
