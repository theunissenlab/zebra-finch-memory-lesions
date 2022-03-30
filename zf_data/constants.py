"""Constants for generating figures in the notebooks.
"""

HVC_COLOR = "#e6438c"  #@param {type: "string"}
CTRL_COLOR = "#777777"  #@param {type: "string"}
NCM_COLOR = "#19b382"  #@param {type: "string"}
NEUTRAL_COLOR = "#1968c2"  #@param {type: "string"}
HVC_LINESTYLE = "--"  #@param {type: "string"}
CTRL_LINESTYLE = (0, (3, 1, 1, 1)) #@param {type: "raw"}
NCM_LINESTYLE = "-"  #@param {type: "string"}
AX_COLOR = "#666666"  #@param {type: "string"}
AXIS_SIZE = 14  #@param {type: "integer"}
LABEL_SIZE = 16  #@param {type: "integer"}

# Used for k <= K_MAX_INITIAL in the paper
K_MAX_INITIAL = 3 #@param {type: "integer"}

# Colors used for coloring by lesion group and linestyles for lesion group
COLORMAP = {
    "NCM": NCM_COLOR,
    "HVC": HVC_COLOR,
    "CTRL": CTRL_COLOR
}

LINEMAP = {
    "NCM": NCM_LINESTYLE,
    "HVC": HVC_LINESTYLE,
    "CTRL": CTRL_LINESTYLE,
}
