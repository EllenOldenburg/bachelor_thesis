import matplotlib
import cycler
import itertools

def set(context="notebook",font_scale=1):
    """ Sets the matplotlib rcParams.
    Available contexts: "paper", "notebook", "talk"
    """

    ##########################################################################
    # Plot sizes
    ##########################################################################
    contexts = ["paper", "notebook", "talk"]

    if context not in contexts:
        raise ValueError("context must be in %s" % ", ".join(contexts))

    # Define colors
    light_gray = ".8"
    dark_gray = ".15"

    # Set up dictionary of default parameters
    base_context = {
        "font.size": 12,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "grid.linewidth": 1,
        "lines.linewidth": 1.75,
        "patch.linewidth": .3,
        "lines.markersize": 7,
        "lines.markeredgewidth": 0,
        "xtick.major.width": 1,
        "ytick.major.width": 1,
        "xtick.minor.width": .5,
        "ytick.minor.width": .5,
        "xtick.major.pad": 7,
        "ytick.major.pad": 7,
        }

    # Scale all the parameters by the same factor depending on the context
    scaling = dict(paper=1.5, notebook=1, talk=1.6,)[context]
    context_dict = {k: v * scaling for k, v in base_context.items()}

    # Set plot size depending on context
    context_dict.update({'figure.figsize': [16,9]})

    # Now independently scale the fonts
    font_keys = ["axes.labelsize", "axes.titlesize", "legend.fontsize",
                 "xtick.labelsize", "ytick.labelsize", "font.size"]
    font_dict = {k: context_dict[k] * font_scale for k in font_keys}
    context_dict.update(font_dict)

    ##########################################################################
    # Style
    ##########################################################################

    # Common parameters
    style_dict = {
        "figure.facecolor": "white",
        "text.color": dark_gray,
        "axes.labelcolor": dark_gray,
        "legend.frameon": False,
        "legend.numpoints": 1,
        "legend.scatterpoints": 1,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.color": dark_gray,
        "ytick.color": dark_gray,
        "axes.axisbelow": True,
        "image.cmap": "rainbow",
        "font.family": ["sans-serif"],
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans",
                            "Bitstream Vera Sans", "sans-serif"],
        "grid.linestyle": "-",
        "lines.solid_capstyle": "round",
        }

    # Add white grid
    style_dict.update({
        "axes.grid": True,
        'axes.facecolor': 'white',
        'axes.edgecolor': 'white',
        "axes.linewidth": 1,
        "grid.color": "#EAEAF2",
        })

    # Ticks
    style_dict.update({
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        })

    ##########################################################################
    # Set color palette and line styles
    ##########################################################################
    # Material colors, based on https://www.materialui.co/colors
    colors = [
        (0.12941176470588237, 0.5882352941176471, 0.9529411764705882),
        (0.0, 0.7372549019607844, 0.8313725490196079),
        (0.0, 0.5882352941176471, 0.5333333333333333),
        (0.2980392156862745, 0.6862745098039216, 0.3137254901960784),
        (0.803921568627451, 0.8627450980392157, 0.2235294117647059),
        (1.0, 0.596078431372549, 0.0),
        (1.0, 0.3411764705882353, 0.13333333333333333),
        (0.9137254901960784, 0.11764705882352941, 0.38823529411764707),
        (0.611764705882353, 0.15294117647058825, 0.6901960784313725)]

    linestyles = ['-', '--', '-.', ':']
    ls = []
    for i in linestyles:
        for j in itertools.repeat(i,len(colors)):
            ls.append(j)

    style_dict.update({
        'axes.prop_cycle': cycler.cycler('color', colors)*len(linestyles) + cycler.cycler("linestyle",ls)
    })

    matplotlib.rcParams.update(context_dict)
    matplotlib.rcParams.update(style_dict)
