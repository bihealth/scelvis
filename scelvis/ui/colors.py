"""Definition of colors."""

tab10 = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]
tab20 = [
    "#1f77b4",
    "#aec7e8",
    "#ff7f0e",
    "#ffbb78",
    "#2ca02c",
    "#98df8a",
    "#d62728",
    "#ff9896",
    "#9467bd",
    "#c5b0d5",
    "#8c564b",
    "#c49c94",
    "#e377c2",
    "#f7b6d2",
    "#7f7f7f",
    "#c7c7c7",
    "#bcbd22",
    "#dbdb8d",
    "#17becf",
    "#9edae5",
]
tab40 = [
    "#393b79",
    "#5254a3",
    "#6b6ecf",
    "#9c9ede",
    "#637939",
    "#8ca252",
    "#b5cf6b",
    "#cedb9c",
    "#8c6d31",
    "#bd9e39",
    "#e7ba52",
    "#e7cb94",
    "#843c39",
    "#ad494a",
    "#d6616b",
    "#e7969c",
    "#7b4173",
    "#a55194",
    "#ce6dbd",
    "#de9ed6",
    "#3182bd",
    "#6baed6",
    "#9ecae1",
    "#c6dbef",
    "#e6550d",
    "#fd8d3c",
    "#fdae6b",
    "#fdd0a2",
    "#31a354",
    "#74c476",
    "#a1d99b",
    "#c7e9c0",
    "#756bb1",
    "#9e9ac8",
    "#bcbddc",
    "#dadaeb",
    "#636363",
    "#969696",
    "#bdbdbd",
    "#d9d9d9",
]

# Set3=['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69',
#      '#fccde5','#d9d9d9','#bc80bd']

my_gradients = [
    [[0, "rgb(180,180,180)"], [1, "rgb(31,119,180)"]],  # grey to blue
    [[0, "rgb(180,180,180)"], [1, "rgb(214,39,40)"]],  # grey to red
    [[0, "rgb(180,180,180)"], [1, "rgb(44,160,44)"]],  # grey to green
    [[0, "rgb(180,180,180)"], [1, "rgb(0,0,0)"]],  # grey to black
]


def get_cm(values):
    """Return color mapping"""
    if len(values) < 10:
        return tab10
    elif len(values) < 20:
        return tab20
    else:
        return tab40
