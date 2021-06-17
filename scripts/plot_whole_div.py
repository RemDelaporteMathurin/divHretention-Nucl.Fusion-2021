'''
Defines a function colorline that draws a (multi-)colored 2D line with coordinates x and y.
The color is taken from optional data in z, and creates a LineCollection.
z can be:
- empty, in which case a default coloring will be used based on the position along the input arrays
- a single number, for a uniform color [this can also be accomplished with the usual plt.plot]
- an array of the length of at least the same length as x, to color according to this data
- an array of a smaller length, in which case the colors are repeated along the curve
The function colorline returns the LineCollection created, which can be modified afterwards.
See also: plt.streamplot
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm
import matplotlib.patches as patches
from divHretention import Exposition

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', size=12)


def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.
    The angle should be given in radians.
    Args:
        origin (float, float): coordinates of origin point
        point (float, float): coordinates of point to be rotated
        angle (float): rotaton angle in radians (counterclockwise)
    Returns:
        float, float: rotated point coordinates
    """
    ox, oy = origin
    px, py = point

    qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
    qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
    return qx, qy


# Data manipulation:

def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments


# Interface to LineCollection:

def colorline(x, y, z=None, cmap=plt.get_cmap('viridis'), norm=plt.Normalize(0.0, 1.0), linewidth=12, alpha=1):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, norm=norm, linewidth=linewidth, alpha=alpha, cmap=cmap)
    ax = plt.gca()
    ax.add_collection(lc)

    return lc


filename = "../data/exposure_conditions_divertor/WEST/West-LSN-P4.4e+21-IP2.000MW.csv"
res = Exposition(filename, filetype="WEST")
res.compute_inventory(1e7)


H = -0.58 - (-0.798)
L = 2.446 - 1.909
angle_div = np.arctan(-H/L)

x = res.arc_length
x = x - x.min()  # starts at 0
x = x*np.cos(angle_div)  # projects on the cartesian plan
x += 1.91  # starts at 1.91

plt.figure(1)
plt.gca().set_aspect('equal')
plt.figure(2)
plt.gca().set_aspect('equal')

thetas = np.linspace(0, 2*np.pi, num=480)

for i, theta in enumerate(thetas):
    points = []
    for point_x in x:
        points.append(rotate(origin=(0, 0), point=(point_x, 0), angle=theta))
    points = np.array(points)
    plt.figure(1)
    lc_inv = colorline(
        points[:, 0], points[:, 1], res.inventory,
        norm=LogNorm(res.inventory.min(), 1e22),
        linewidth=4)
    plt.figure(2)
    lc_T = colorline(
        points[:, 0], points[:, 1], res.temperature,
        norm=plt.Normalize(res.temperature.min(), res.temperature.max()),
        cmap=plt.get_cmap('inferno'),
        linewidth=4)

plt.figure(1)
cb = plt.colorbar(lc_inv)
cb.set_label("Inventory (H/m)")

plt.figure(2)
cb = plt.colorbar(lc_T)
cb.set_label("$T$ (K)")

plt.figure(1)
xmax = 2.5
plt.xlim(-xmax, xmax)
plt.ylim(-xmax, xmax)

plt.xlabel("x (m)")
plt.ylabel("y (m)")

# plt.axis('off')

plt.figure(2)
xmax = 2.5
plt.xlim(-xmax, xmax)
plt.ylim(-xmax, xmax)

plt.xlabel("x (m)")
plt.ylabel("y (m)")

# plt.axis('off')
plt.figure(1)
plt.savefig("../Figures/WEST/circle_inventory.pdf")
plt.show()
