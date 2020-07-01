import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


class color_by_reward(object):
    @staticmethod
    def get(x):
        if "Rewarded" in x:
            return "#0AA5D8"
        else:
            return "#C62533"


def plot_data(
        df,
        grouping,
        force_len=None,
        index_by="time",
        label_order=None,
        label_to_color=None,
        tick_height=0.1,
        figsize=None
    ):

    # Reindex the dataframe for plotting
    original_index = df.index
    df.index = pd.Series(np.arange(len(df)))
    groupings = list(df.groupby(grouping))

    n_categories = len(groupings)

    if figsize is None:
        fig = plt.figure(
            facecolor="white",
            edgecolor="white",
            figsize=(10, 4 + 5 * tick_height * n_categories)
        )
    else:
        fig = plt.figure(
            facecolor="white",
            edgecolor="white",
            figsize=figsize
        )

    events_ax = fig.gca()
    prob_ax = events_ax.twinx()

    events_ax.set_ylim(-0.2 - tick_height * n_categories, 1.2 + tick_height * n_categories)
    prob_ax.set_ylim(-0.2 - tick_height * n_categories, 1.2 + tick_height * n_categories)

    for group_idx, (group_keys, group_df) in enumerate(groupings):
        # Interrupted trials will be plotted on top of plot (starting at y=1)
        # going in the positive direction
        # Non-interrupted trials will be plotted on bottom of plot (starting at y=0)
        # going in the negative direction
        interrupted = group_df["Interrupt"].apply(lambda x: 1 if x else 0)
        increment_direction = group_df["Interrupt"].apply(lambda x: 1 if x else -1)

        # Plot event tick marks
        scatter_plot = events_ax.scatter(
            group_df.index,
            (
                ((1 * interrupted) + (2 * tick_height * increment_direction)) +  # tick base position
                increment_direction * tick_height * group_idx                    # offset each group
            ),
            s=50,
            marker="|",
            color=color_by_reward.get(group_keys),
            label=" ".join(group_keys)
        )

        events_ax.vlines(x=0, ymin=1 + tick_height, ymax=1.3 + tick_height * n_categories, color='black', linewidth=2)
        events_ax.vlines(x=0, ymin=-0.3 - tick_height * n_categories, ymax=-tick_height, color='black', linewidth=2)

        # Plot a line showing windowed probability of interruption
        win_size=20
        win_size_half = win_size // 2
        rolled = group_df["Interrupt"].rolling(win_size, center=True).mean()

        # Fill in nans at beginning/end by the first/last value
        if len(rolled) > win_size_half:
            rolled.iloc[:win_size_half] = rolled.iloc[win_size_half]
            rolled.iloc[-win_size_half:] = rolled.iloc[-win_size_half - 1]

        prob_ax.plot(
            group_df.index,
            rolled,
            label=scatter_plot.get_label(),
            alpha=1.0,
            linewidth=3,
            color=scatter_plot.get_edgecolor()[0]
        )


    # Draw borders between probability plot and trial ticks, stylize by shading background
    events_ax.hlines([-0.01, 1.01], *events_ax.get_xlim(), linewidth=2, linestyle=":", color="Grey")
    events_ax.fill_between(events_ax.get_xlim(), [-0.01, -0.01], [1.01, 1.01], color="0.95", zorder=0)

    prob_ax.set_xlim(0, force_len or len(df))

    # Clean up and label axes
    events_ax.xaxis.set_tick_params(labelsize=16)
    events_ax.set_yticks([0, 1])
    events_ax.set_yticklabels([0.0, 1.0], size=16)
    events_ax.set_ylabel("Prob.\ninterrupt", fontsize=16)
    events_ax.set_xticks([])
    events_ax.set_xlabel("Trial", fontsize=16)
    events_ax.set_yticks([0.2, 0.4, 0.6, 0.8], minor=True)
    events_ax.grid(which='minor', alpha=0.8, linestyle=":")
    prob_ax.set_yticks([])

    events_ax.spines['top'].set_visible(False)
    events_ax.spines['right'].set_visible(False)
    events_ax.spines['bottom'].set_visible(False)
    events_ax.spines['left'].set_visible(False)
    prob_ax.spines['top'].set_visible(False)
    prob_ax.spines['right'].set_visible(False)
    prob_ax.spines['bottom'].set_visible(False)
    prob_ax.spines['left'].set_visible(False)

    # Label tick marks
    events_ax.text(0, 1 + (2 * tick_height) + tick_height * 0.5 * n_categories, "Int.  ", fontsize=16, horizontalalignment="right", verticalalignment="center")
    events_ax.text(0, -(2 * tick_height) - tick_height * 0.5 * n_categories, "Wait  ",  fontsize=16, horizontalalignment="right", verticalalignment="center")

    events_ax.vlines(x=0, ymin=0, ymax=1, color='black', linewidth=2)

    df.index = original_index

    return fig


def set_oddsratio_yticks(ax, biggest, smallest=None):
    """Determine and set the yticks of an axis given the data range

    Generates a pleasant set of ytick labels and spacing for a given
    range of odds ratios.

    Typecasts the y values into multiples (e.g. x1, x2, x4, etc) when the
    odds ratio is > 1 or as fractions when the odds ratio is less than 1
    (e.g. x1/2, x1/4, etc).

    It ensures that:
    * only powers of 2 are shown
    * y=1 is always labeled
    * there is a maximum of 5 y-values labeld
    """
    if not smallest:
        smallest = -biggest
    if smallest >= -1:
        smallest = -1

    ax.set_yscale("log")

    abs_biggest = max(np.abs(smallest), np.abs(biggest))

    powers = np.arange(0, abs_biggest + 1)
    n = len(powers)
    powers = powers[::n // 6 + 1]
    vals = np.concatenate([-powers, powers[1:]])

    ticks = np.power(2., vals)
    labels = [r"x{:d}".format(int(2 ** v)) if v >= 0 else r"x1/{:d}".format(int(2 ** -v)) for v in vals]

    ax.set_ylabel("Odds Ratio", fontsize=16)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=16)
    ax.hlines(1, *plt.xlim(), linestyle="--", zorder=-1)
    ax.set_ylim(np.power(2., smallest), np.power(2., biggest))
