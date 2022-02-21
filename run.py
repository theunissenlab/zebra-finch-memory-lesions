import matplotlib.pyplot as plt
import pandas as pd

from zf_data import Tsvk


if __name__ == "__main__":
    df = pd.read_csv("~/Desktop/TrialDataNew.csv")
    base = Tsvk(df[(df.VocalizerSet == "S1") & (df.LesionStage == "prelesion")], k_max=21)
    result = base.logOR()

    plt.fill_between(result["k"], result["logOR"] - 2 * result["SE"], result["logOR"] + 2 * result["SE"], alpha=0.5)
    plt.plot(result["logOR"], label="Prelesion")

    for treatment in ["NCM", "HVC", "CTRL"]:
        T = Tsvk(df[(df.VocalizerSet == "S1") & (df.LesionStage == "postlesion") & (df.SubjectTreatment == treatment)], k_max=11)
        result = T.logOR()
        plt.plot(result["logOR"], label=treatment)

    plt.ylim(-1, 8)

    plt.legend()
    plt.show()
