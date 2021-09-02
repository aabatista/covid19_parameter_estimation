import numpy as np
import matplotlib.collections as mcoll
from matplotlib import cm
def colorline(ax, x, y, z, cmap=cm.get_cmap('copper'), linewidth=3, alpha=1.0):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Specify colors in the array z
    """
    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, linewidth=linewidth
        , alpha=alpha)
    ax.add_collection(lc)
    return lc

def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    # points = np.array([x, y], dtype=float).T.reshape(-1, 1, 2)
    #segments = np.concatenate([points[:-1], points[1:]], axis=1)
    coords = list(zip(x, y))
    protoSegs = list(zip(coords[:-1], coords[1:]))
    segments = [(start, end) for start, end in protoSegs]
    return segments

