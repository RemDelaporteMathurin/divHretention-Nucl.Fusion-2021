import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import numpy as np
from scipy.stats import linregress

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from divHretention import data as data_module
from divHretention import reflection_coeff

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)


with pkg_resources.path(data_module, "data_TRIM_energy_angle.csv") as p:
    data = np.genfromtxt(p, delimiter=";", names=True)


plt.figure(1, figsize=[6.4, 4])
angles = np.unique(data["theta_inc"])
colormap = cm.Blues
for angle in angles:
    indexes = np.where(data["theta_inc"] == angle)

    energy = data["Incident_energy"][indexes]
    r = data["Reflection_coeff"][indexes]

    plt.figure(1)
    offset = 50
    colour = colormap((angle + offset)/(90 + offset))
    plt.scatter(
        energy,
        r,
        label="{:.0f} °".format(angle),
        alpha=0.6,
        color=colour)
    plt.plot(
        energy,
        [reflection_coeff(energy[i], angle) for i in range(len(energy))],
        color=colour)
    plt.annotate("{:.0f} °".format(angle), (energy[-1] * 1.05, r[-1]), color=colour)
    plt.figure(2)

    stdev_x = .005 * (max(angles) - min(angles))
    plt.scatter(
        np.ones(energy.shape)*angle + np.random.randn(len(energy)) * stdev_x,
        r,
        c=energy,
        cmap="viridis",
        norm=LogNorm(),
        alpha=0.6)

plt.figure(1)
plt.ylim(bottom=0, top=1)
plt.xscale("log")
plt.xlabel("Incident energy (eV)")
plt.ylabel("Reflection coefficient")
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

plt.savefig("../Figures/reflection_coeff.pdf")
plt.savefig("../Figures/reflection_coeff.svg")

plt.figure(2)
plt.colorbar()
plt.ylim(bottom=0, top=1)
plt.xlabel("Angle of incidence (°)")
plt.ylabel("Reflection coefficient")
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

plt.show()
