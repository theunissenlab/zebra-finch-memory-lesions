import os

import numpy as np
import pandas as pd


if os.environ.get("DATADIR"):
    DATADIR = os.environ["DATADIR"]
else:
    _CODEDIR, _ = os.path.split(__file__)
    DATADIR = os.path.join(_CODEDIR, "..", "data")


def load_trials(limit_rows=None):
    """Load pandas DataFrame from TrialData.csv"""
    return pd.read_csv(
        os.path.join(DATADIR, "behavior", "TrialData.csv"),
        parse_dates=["Time"],
        converters={
            "Date": lambda d: pd.to_datetime(d).date(),
            "RT": pd.to_timedelta
        },
        nrows=limit_rows
    )



