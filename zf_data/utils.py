import matplotlib as mpl

from .constants import AX_COLOR


def parse_p(p_value: float):
    """Round p-value to 3 digits, or return "< 0.001" for smaller values
    
    Edit this function if submission requirements need different p-value formatting.
    """
    if p_value < 0.001:
        sig_str = " (***)"
    elif p_value < 0.01:
        sig_str = " (**)"
    elif p_value < 0.05:
        sig_str = " (*)"
    else:
        sig_str = ""
        
    if p_value < 0.0005:
        return f"p < 0.001{sig_str}"
    else:
        return f"p = {p_value:.3f}{sig_str}"


def setup_mpl_params():
    """Set all the mpl axes colors to a default greyish color
    """
    mpl.rcParams["grid.color"] = AX_COLOR
    mpl.rcParams["axes.edgecolor"] = AX_COLOR
    mpl.rcParams["xtick.labelcolor"] = AX_COLOR
    mpl.rcParams["ytick.labelcolor"] = AX_COLOR
    mpl.rcParams["xtick.color"] = AX_COLOR
    mpl.rcParams["ytick.color"] = AX_COLOR
    mpl.rcParams["ytick.color"] = AX_COLOR

    mpl.rcParams["axes.titlecolor"] = AX_COLOR
    mpl.rcParams["axes.labelcolor"] = AX_COLOR
    mpl.rcParams["figure.edgecolor"] = AX_COLOR
    mpl.rcParams["grid.color"] = AX_COLOR
    mpl.rcParams["legend.labelcolor"] = AX_COLOR
    mpl.rcParams["legend.edgecolor"] = AX_COLOR
    mpl.rcParams["text.color"] = AX_COLOR

