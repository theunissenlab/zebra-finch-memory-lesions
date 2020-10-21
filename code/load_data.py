import os

import numpy as np
import pandas as pd
import wavio


CODEDIR, _ = os.path.split(__file__)
DATADIR = os.path.join(CODEDIR, "..")


def load_behavioral_data():
    """Load pandas DataFrame from TrialData.csv"""
    return pd.read_csv(
        os.path.join(DATADIR, "behavior", "TrialData.csv"),
        parse_dates=["Time"],
        converters={
            "Date": lambda d: pd.to_datetime(d).date(),
            "RT": pd.to_timedelta
        }
    )


def load_lesion_data():
    """Load pandas DataFrame from TrialData.csv"""
    return pd.read_csv(
        os.path.join(DATADIR, "behavior", "LesionData.csv"),
    )


def load_spike_data():
    df = pd.read_pickle(os.path.join("data", "ephys", "UnitData.pkl"))
    df["n_trials"] = df.apply(lambda row: np.sum(row["good_trials"]), axis=1)
    df["spike_times"] = df.apply(lambda row: np.array(row["spike_times"])[np.array(row["good_trials"])], axis=1)
    df = df.drop("good_trials", axis=1)
    df = df.set_index(["unit_id", "file"]).sort_index()
    return df


def load_stim_data():
    df = pd.read_pickle(os.path.join("data", "ephys", "StimData.pkl")
    return df


def load_behavior_stimulus(stim_file):
    """Return wavio.Wav object for given stimulus filename"""
    return wavio.read(os.path.join("data", "behavior", "stimuli", stim_file))


def load_ephys_stimulus(stim_file):
    """Return wavio.Wav object for given stimulus filename"""
    return wavio.read(os.path.join("data", "ephys", "stimuli", stim_file))


def load_unit_waveforms(unit_id):
    return np.load(os.path.join("data", "ephys", "unit_waveforms", unit_id))


__all__ = [
    "load_behavioral_data",
    "load_lesion_data",
    "load_spike_data",
    "load_stim_data",
    "load_behavior_stimulus",
    "load_ephys_stimulus",
    "load_unit_waveforms"
]
