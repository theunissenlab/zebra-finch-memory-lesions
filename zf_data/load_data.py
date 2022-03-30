import datetime
import glob
import os

import numpy as np
import pandas as pd


if os.environ.get("DATADIR"):
    DATADIR = os.environ["DATADIR"]
else:
    _CODEDIR, _ = os.path.split(__file__)
    DATADIR = os.path.join(_CODEDIR, "..", "data")


EXCLUSION_SUBJECTS = [
    "BluWhi3230M",
    "GreWhi2703M",
]
EXCLUSION_DATE = datetime.date(2020, 11, 20)
EXCLUSION_TIME = datetime.datetime(2020, 11, 20, 12, 0, 0)


def load_trials(valid_only=True, limit_rows=None):
    """Load pandas DataFrame from TrialData.csv"""
    df = pd.read_csv(
        os.path.join(DATADIR, "behavior", "TrialData.csv"),
        parse_dates=["Time"],
        converters={
            "Date": lambda d: pd.to_datetime(d).date(),
            "RT": pd.to_timedelta
        },
        nrows=limit_rows
    )

    if valid_only:
        # Filter out data where audio output was not working
        # (2 subjects on 2020-11-20 after 12pm)
        return df[
            ~(
                df.Subject.isin(EXCLUSION_SUBJECTS) 
                & (df.Date == EXCLUSION_DATE)
                & (df.Time > EXCLUSION_TIME)
            )
        ]
    else:
        return df


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
