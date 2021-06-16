import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from divHretention import data as data_module

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)


with pkg_resources.path(data_module, "data_TRIM_energy_angle.csv") as p:
    data = np.genfromtxt(p, delimiter=";", names=True)


plt.figure(figsize=[6.4, 4])

for angle in np.unique(data["theta_inc"]):
    indexes = np.where(data["theta_inc"] == angle)

    energy = data["Incident_energy"][indexes]
    R_p = data["Implantation_range"][indexes] * 1e9

    # add jitter for better visualisation
    stdev_x = .0045 * (max(np.log10(energy)) - min(np.log10(energy)))
    stdev_y = .0045 * (max(R_p) - min(R_p))
    plt.scatter(
        10**(np.log10(energy) + np.random.randn(len(energy)) * stdev_x),
        R_p + np.random.randn(len(R_p)) * stdev_y,
        label="{:.0f} °".format(angle),
        alpha=0.6)

# add regression
energy_all = data["Incident_energy"]
R_p_all = data["Implantation_range"] * 1e9
res = linregress(np.log10(energy_all), np.log10(R_p_all))
a = 10**res.intercept
b = res.slope

x = np.logspace(0, np.log10(2e3))
plt.plot(x, a*x**b, linestyle="dashed")

annotation = "$" + "{:.2f}".format(a) + "x ^ {" + "{:.2f}".format(b) + "}$"
plt.annotate(annotation, (x[-1], a*x[-1]**b), color="tab:blue")

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

plt.savefig("../Figures/implantation_range.pdf")
plt.savefig("../Figures/implantation_range.svg")

plt.show()
