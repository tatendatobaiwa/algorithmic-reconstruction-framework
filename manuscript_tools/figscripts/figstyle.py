"""Shared figure style for the manuscript figures (MDPI: >= 600 dpi, RGB, PNG).
Okabe-Ito colour-blind-safe palette; serif text to match the manuscript."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OKABE = {"orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
         "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7", "grey": "#7F7F7F", "black": "#000000"}
SELECTOR = {"ECC": OKABE["vermillion"], "Coverage-first": OKABE["blue"], "Time-aware": OKABE["green"],
            "FCFS": OKABE["purple"], "Exact": OKABE["black"]}
SEASON = {"dry": OKABE["orange"], "wet": OKABE["blue"]}
FULL_W, HALF_W = 7.0, 3.4     # inches
DPI = 600


def setup():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Liberation Serif", "DejaVu Serif", "Times New Roman"],
        "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5, "legend.fontsize": 7.5,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "savefig.dpi": DPI, "savefig.bbox": "tight", "savefig.pad_inches": 0.03, "figure.dpi": 100,
        "axes.titleweight": "bold", "axes.titlelocation": "left"})


def save(fig, path):
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    from PIL import Image
    im = Image.open(path)
    assert im.mode in ("RGB", "RGBA") and abs(im.info.get("dpi", (0, 0))[0] - DPI) < 1, (im.mode, im.info.get("dpi"))
    return im.size
