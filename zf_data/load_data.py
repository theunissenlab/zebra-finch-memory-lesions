import glob
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


def load_lesion_summary_table():
    """Load pandas DataFrame from TrialData.csv"""
    return pd.read_csv(
        os.path.join(DATADIR, "behavior", "LesionQuantificationSummary.csv"),
    )


def load_subject_lesion_tables():
    tables = {}
    for csv_file in glob.glob(os.path.join(DATADIR, "behavior", "LesionQuantification", "*.csv")):
        subject = os.path.splitext(os.path.basename(csv_file))[0]
        tables[subject] = pd.read_csv(csv_file)

    return tables
