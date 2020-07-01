import os

import pandas as pd
import wavio


CODEDIR, _ = os.path.split(__file__)
DATADIR = os.path.join(CODEDIR, "..")


def load_data():
    """Load pandas DataFrame from TrialData.csv"""
    return pd.read_csv(
        os.path.join(DATADIR, "TrialData.csv"),
        parse_dates=["Time"],
        converters={
            "Date": lambda d: pd.to_datetime(d).date(),
            "RT": pd.to_timedelta
        }
    )


def load_stimulus(stim_file):
    """Return wavio.Wav object for given stimulus filename"""
    return wavio.read(os.path.join("..", "stimuli", stim_file))


__all__ = ["load_data", "load_stimulus"]
