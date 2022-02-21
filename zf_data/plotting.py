from itertools import product

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


def border(ax=None, left=False, right=False, top=False, bottom=False):
    if ax is None:
        ax = plt.gca()
    ax.spines["left"].set_visible(left)
    ax.spines["right"].set_visible(right)
    ax.spines["top"].set_visible(top)
    ax.spines["bottom"].set_visible(bottom)


def _get_layout(n, max_columns=None):
    w = 1
    h = 1

    increment_width = False
    while w * h < n:
        if increment_width and (not max_columns or w < max_columns):
            w += 1
        else:
            h += 1
        increment_width = not increment_width

    return w, h


def fig_grid(n, ax_size=(3, 2), max_columns=None, columns=None, xpad=0.0, ypad=0.0, start="top"):
    """Create a grid of n axes for plotting

    Tries its best to keep the number of rows and columns even, up to a max number of columns

    Params
    ======
    n: number of subplots to create
    ax_size (tuple): width and height of each axis to create
    max_columns: tries to keep number of rows and columns even, but will not go
        above max_columns
    xpad (float, default=0): How much padding between figure columns
    ypad (float, default=0): How much padding between figure rows
    start (str, either "top" or "bottom"): Iterate over axes starting from top or from bottom

    Returns
    =======
    Returns a Figure refernce and a list of n matplotlib Axes objects
    """
    # dimension
    if columns is not None:
        w = columns
        h = n // columns + 1 * bool(n % columns)
    else:
        w, h = _get_layout(n, max_columns=max_columns)

    fig = plt.figure(figsize=ax_size)
    axes = []

    # Start from top
    col_positions = [i * (1 + xpad) for i in range(w)]
    row_positions = [i * (1 + ypad) for i in range(h)]
    if start == "top":
        row_positions = row_positions[::-1]
    elif start == "bottom":
        pass

    ax_positions = list(product(col_positions, row_positions))[:n]
    for pos_x, pos_y in ax_positions:
        axes.append(fig.add_axes([pos_x + 0.1, pos_y + 0.1, 0.8, 0.8]))

    return fig, axes
    curr_x = 0
    curr_y = 0
    for i in range(n):
        axes.append(fig.add_axes([curr_x + 0.1, curr_y + 0.1, 0.8, 0.8]))
        curr_x += 1
        if curr_x == w:
            curr_x = 0
            curr_y += 1

    return fig, axes


def plot_pecking_test_data(
        df,
        grouping,
        force_len=None,
        index_by="time",
        label_order=None,
        label_to_color=None,
        tick_height=0.1,
        ticks=True,
        reward_to_color=None,
        figsize=None,
        linekwargs=None,
        mark_days=False,
    ):
    """Plot behavioral data from pecking test
    """
    if reward_to_color is None:
        reward_to_color = color_by_reward
    if linekwargs is None:
        linekwargs = {}

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
    
    if ticks:
        events_ax.set_ylim(-0.2 - tick_height * n_categories, 1.2 + tick_height * n_categories)
        prob_ax.set_ylim(-0.2 - tick_height * n_categories, 1.2 + tick_height * n_categories)
    else:
        events_ax.set_ylim(-0.0, 1)
        prob_ax.set_ylim(-0.0, 1)
        
    for group_idx, (group_keys, group_df) in enumerate(groupings):
        # Interrupted trials will be plotted on top of plot (starting at y=1)
        # going in the positive direction
        # Non-interrupted trials will be plotted on bottom of plot (starting at y=0)
        # going in the negative direction
        
        interrupted = group_df["Interrupt"].apply(lambda x: 1 if x else 0)
        increment_direction = group_df["Interrupt"].apply(lambda x: 1 if x else -1)

        # Plot event tick marks
        if ticks:
            scatter_plot = events_ax.scatter(
                group_df.index,
                (
                    ((1 * interrupted) + (2 * tick_height * increment_direction)) +  # tick base position
                    increment_direction * tick_height * group_idx                    # offset each group
                ),
                s=50,
                marker="|",
                color=reward_to_color.get(group_keys),
                label=" ".join(group_keys)
            )

            events_ax.vlines(x=0, ymin=1 + tick_height, ymax=1.3 + tick_height * n_categories, color='black', linewidth=2)
            events_ax.vlines(x=0, ymin=-0.3 - tick_height * n_categories, ymax=-tick_height, color='black', linewidth=2)

        # Plot a line showing windowed probability of interruption
        win_size=20 if group_keys[0] == "Nonrewarded" else 12
        win_size_half = win_size // 2
        rolled = group_df["Interrupt"].rolling(win_size, center=True).mean()

        # Fill in nans at beginning/end by the first/last value
        if len(rolled) > win_size_half:
            rolled.iloc[:win_size_half] = rolled.iloc[win_size_half]
            rolled.iloc[-win_size_half:] = rolled.iloc[-win_size_half - 1]

        kwargs = dict(
            label=scatter_plot.get_label() if ticks else None,
            alpha=1.0,
            linewidth=3,
            color=reward_to_color.get(group_keys),
        )
        kwargs.update(linekwargs)
        
        prob_ax.plot(
            group_df.index,
            rolled,
            **kwargs
        )
        
    if mark_days:
        import datetime
        # Draw vertical lines where trials transition days
        idx = np.where(df["Date"].diff() >= datetime.timedelta(days=1))[0]
        prob_ax.vlines(idx, *prob_ax.get_ylim(), linestyle="-", color="Black", zorder=100, linewidth=0.5)


    # Draw borders between probability plot and trial ticks, stylize by shading background
    if ticks:
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
    if ticks:
        events_ax.text(0, 1 + (2 * tick_height) + tick_height * 0.5 * n_categories, "Int.  ", fontsize=16, horizontalalignment="right", verticalalignment="center")
        events_ax.text(0, -(2 * tick_height) - tick_height * 0.5 * n_categories, "Wait  ",  fontsize=16, horizontalalignment="right", verticalalignment="center")

        events_ax.vlines(x=0, ymin=0, ymax=1, color='black', linewidth=2)

    df.index = original_index

    return fig


def set_oddsratio_yticks(ax, biggest, smallest=None, convert_log=True, set_label=True):
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

    if convert_log:
        ax.set_yscale("log")

    abs_biggest = max(np.abs(smallest), np.abs(biggest))

    powers = np.arange(0, abs_biggest + 1)
    n = len(powers)
    powers = powers[::n // 6 + 1]
    vals = np.concatenate([-powers, powers[1:]])
    vals = vals[(vals >= smallest) & (vals <= biggest)]

    if convert_log:
        ticks = np.power(2., vals)
    else:
        ticks = vals
        
    labels = [r"x{:d}".format(int(2 ** v)) if v >= 0 else r"x1/{:d}".format(int(2 ** -v)) for v in vals]

    if set_label:
        ax.set_ylabel("Odds Ratio", fontsize=16)

    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=16)

    if convert_log:
        ax.hlines(1, *plt.xlim(), linestyle="--", zorder=-1)
        ax.set_ylim(np.power(2., smallest), np.power(2., biggest))
