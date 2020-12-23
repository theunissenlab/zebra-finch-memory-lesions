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
    lesion_data = pd.read_csv(
        os.path.join(DATADIR, "behavior", "LesionData.csv"),
        converters={
            "Date_Lesioned": lambda d: pd.to_datetime(d).date()
        }
    )

    if len(lesion_data) != 10:
        raise Exception("The lesion data csv has been updated. You must update this function to drop"
                       " the appropriate rows.")

    lesion_data = lesion_data.drop(index=6).reindex()
    return lesion_data


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


def load_auditory_info():
    if not os.path.exists(os.path.join(DATADIR, "ephys", "AuditoryUnitTable.pkl")):
        raise IOError("AuditoryUnitTable.pkl does not exist. May need to run notebook EPHYS3. Identify Auditory Units")

    df = pd.read_pickle(os.path.join(DATADIR, "ephys", "AuditoryUnitTable.pkl"))
    return df


def load_rendition_data():
    """Loads ephys data split into individual renditions

    Derived from data loaded in load_ephys_data() but seperates each
    sub-rendition of a trial into its own row, providing the columns
    of "rendition_idx", "trial_onset", and "preceding_silence_duration"
    as references to the original trial information.
    """
    df = pd.read_pickle(os.path.join(DATADIR, "ephys", "UnitRenditionData.pkl"))
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
    with np.load(os.path.join(DATADIR, "ephys", "unit_waveforms", "{}.npz".format(unit_id))) as data:
        return {
            "spike_times": data["spike_times"],
            "spike_waveforms": data["spike_waveforms"],
        }


__all__ = [
    "load_behavioral_data",
    "load_lesion_data",
    "load_spike_data",
    "load_rendition_data",
    "load_stim_data",
    "load_behavior_stimulus",
    "load_ephys_stimulus",
    "load_unit_waveforms"
]
