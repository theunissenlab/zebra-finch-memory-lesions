import datetime
import os

import numpy as np
import pandas as pd
import wavio


CODEDIR, _ = os.path.split(__file__)
DATADIR = os.path.join(CODEDIR, "..", "data")


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
        converters={
            "Date_Lesioned": lambda d: pd.to_datetime(d).date()
        }
    )


def load_spike_data():
    df = pd.read_pickle(os.path.join(DATADIR, "ephys", "UnitData.pkl"))

    # TODO eventually apply these to the original dataframe so we dont have to compute it every time we load
    df["n_trials"] = df.apply(lambda row: np.sum(row["good_trials"]), axis=1)
    df["spike_times"] = df.apply(lambda row: np.array(row["spike_times"])[np.array(row["good_trials"])], axis=1)
    df["subject"] = df["site"].apply(lambda site: site.split("_")[0])
    def _parse_date(site):
        date_str = site.split("_")[1]
        year = 2000 + int(date_str[:2])
        month = int(date_str[2:4])
        day = int(date_str[4:])
        return datetime.date(year, month, day)

    df["date"] = df["site"].apply(_parse_date)
    df = df.drop("good_trials", axis=1)
    df = df.sort_values(["unit_id", "stimulus_file"])
    return df


def load_stim_data():
    df = pd.read_pickle(os.path.join(DATADIR, "ephys", "StimData.pkl"))
    return df


def load_behavior_stimulus(stim_file):
    """Return wavio.Wav object for given stimulus filename"""
    return wavio.read(os.path.join(DATADIR, "behavior", "stimuli", stim_file))


def load_ephys_stimulus(stim_file):
    """Return wavio.Wav object for given stimulus filename"""
    return wavio.read(os.path.join(DATADIR, "ephys", "stimuli", stim_file))


def load_unit_waveforms(unit_id):
    return np.load(os.path.join(DATADIR, "ephys", "unit_waveforms", "{}.npz".format(unit_id)))


__all__ = [
    "load_behavioral_data",
    "load_lesion_data",
    "load_spike_data",
    "load_stim_data",
    "load_behavior_stimulus",
    "load_ephys_stimulus",
    "load_unit_waveforms"
]
