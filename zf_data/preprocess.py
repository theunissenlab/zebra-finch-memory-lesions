import pandas as pd


def get_experiment_dates_from_trial_df(df):
    """Generates the ExperimentDates.csv table from the TrialData.csv
    """
    table = dict((subject, {"Subject": subject}) for subject in df.Subject.unique())
    for (subject, voc_set, lesion, ladder), subdf in df.sort_values("Date").groupby(["Subject", "VocalizerSet", "LesionStage", "LadderStage"]):
        # subdf = subdf.sort_values("Date")
        table[subject][f"{voc_set},{lesion},{ladder}"] = "{},{}".format(
            subdf.iloc[0]["Date"].strftime("%Y-%m-%d"),
            subdf.iloc[-1]["Date"].strftime("%Y-%m-%d")
        )

    index = df.Subject.unique()
    result = [table[subject] for subject in index]
    output = pd.DataFrame(result, index=index)

    # Sort them by experiment order
    output = output[[
        "S1,prelesion,DCvsDC_1v1",
        "S1,prelesion,DCvsDC_4v4",
        "S1,prelesion,DCvsDC_6v6_d1",
        "S1,prelesion,DCvsDC_6v6_d2",
        "S1,prelesion,SovsSo_1v1",
        "S1,prelesion,SovsSo_4v4",
        "S1,prelesion,SovsSo_8v8_d1",
        "S1,prelesion,SovsSo_8v8_d2",

        "S1,postlesion,SovsSo_1v1",
        "S1,postlesion,SovsSo_8v8_d2",
        "S1,postlesion,DCvsDC_1v1",
        "S1,postlesion,DCvsDC_6v6_d2",

        "S2,postlesion,DCvsDC_1v1_S2",
        "S2,postlesion,DCvsDC_4v4_S2",
        "S2,postlesion,DCvsDC_6v6_d1_S2",
        "S2,postlesion,DCvsDC_6v6_d2_S2",
        "S2,postlesion,SovsSo_1v1_S2",
        "S2,postlesion,SovsSo_4v4_S2",
        "S2,postlesion,SovsSo_8v8_d1_S2",
        "S2,postlesion,SovsSo_8v8_d2_S2",
    ]]
    output.index.name = "Subject"

    return output

