import matplotlib.pyplot as plt
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

plt.figure(figsize=[6.4, 4])

for angle in np.unique(data["theta_inc"]):
    indexes = np.where(data["theta_inc"] == angle)

    energy = data["Incident_energy"][indexes]
    R_p = data["Reflection_coeff"][indexes]

    # add jitter for better visualisation
    plt.scatter(
        energy,
        R_p,
        label="{:.0f} °".format(angle),
        alpha=0.6)
    plt.plot(
        energy,
        [reflection_coeff(energy[i], angle) for i in range(len(energy))])


plt.ylim(bottom=0)
plt.xscale("log")
# plt.yscale("log")
plt.legend()
plt.xlabel("Incident energy (eV)")
plt.ylabel("Implantation range (nm)")
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

plt.savefig("../Figures/reflection_coeff.pdf")
plt.savefig("../Figures/reflection_coeff.svg")

plt.show()
