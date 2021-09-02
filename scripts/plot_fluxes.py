import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

from divHretention import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)

quantities = [
            "heat_flux",
            "atom_flux", "ion_flux",
            "ion_energy", "atom_energy",
            ]

folder = "../data/exposure_conditions_divertor/WEST/"

# power scan

input_powers = [
    0.449,
    1,
    1.5,
    2
]
puff_rate = 2.5e21
filenames = [
    folder + "West-LSN-P{:.1e}-IP{:.3f}MW.csv".format(puff_rate, input_power)
    for input_power in input_powers]


colormap = cm.inferno
sm = plt.cm.ScalarMappable(cmap=colormap, norm=Normalize(vmin=min(input_powers), vmax=max(input_powers)))
colours = [colormap((IP - min(input_powers))/max(input_powers)) for IP in input_powers]

my_plot = plot_along_divertor(
    filenames=filenames,
    filetypes="WEST",
    quantities=quantities,
    figsize=(6, 5),
    colors=colours)
for ax in my_plot.axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.colorbar(
    sm, label="Input power (MW)",
    ax=my_plot.axs)

custom_lines = [Line2D([0], [0], color="black", linestyle="solid"),
                Line2D([0], [0], color="black", linestyle="dashed")]
my_plot.axs[1].legend(custom_lines, ['Ions', 'Neutrals'], loc="lower right")
my_plot.axs[-1].get_legend().remove()

my_plot.axs[0].set_ylabel("Heat flux \n (W m$^{-2}$)")
my_plot.axs[1].set_ylabel("Incident flux \n (m$^{-2}$ s$^{-1}$)")
my_plot.axs[2].set_ylabel("Incident energy (eV)")

# plt.show()

# density scan
input_power = 0.449
Ps = [
    4.4e20,
    1.0e21,
    1.3e21,
    1.6e21,
    2.0e21,
    2.5e21,
    2.9e21,
    3.3e21,
    3.83e21,
    4.36e21,
    4.7e21,
]
colormap = cm.viridis
sm = plt.cm.ScalarMappable(cmap=colormap, norm=Normalize(vmin=min(Ps), vmax=max(Ps)))
colours = [colormap((P - min(Ps))/max(Ps)) for P in Ps]

filenames = ["../data/exposure_conditions_divertor/WEST/West-LSN-P{:.1e}-IP{:.3}MW.csv".format(P, input_power) for P in Ps]

my_plot = plot_along_divertor(
    filenames=filenames,
    filetypes="WEST",
    quantities=quantities,
    figsize=(6, 5),
    colors=colours)
for ax in my_plot.axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.colorbar(
    sm, label="Puffing rate (molecule s$^{-1}$)",
    ax=my_plot.axs)

custom_lines = [Line2D([0], [0], color="black", linestyle="solid"),
                Line2D([0], [0], color="black", linestyle="dashed")]
my_plot.axs[1].legend(custom_lines, ['Ions', 'Neutrals'], loc="lower right")
my_plot.axs[-1].get_legend().remove()

my_plot.axs[0].set_ylabel("Heat flux \n (W m$^{-2}$)")
my_plot.axs[1].set_ylabel("Incident flux \n (m$^{-2}$ s$^{-1}$)")
my_plot.axs[2].set_ylabel("Incident energy (eV)")

# plt.show()

# ITER case
numbers = [
    2404,
    2403,
    2401,
    2402,
    2399,
    2400,
    2398,
    2397,
    2396,
]

divertor_pressure = [
    11.19639589,
    9.248295137,
    6.889631634,
    8.169794716,
    3.803158809,
    5.132170779,
    2.832874099,
    2.250856857,
    1.752557796,
]
divertor_pressure = np.array(divertor_pressure)
numbers = np.array(numbers)
arr1inds = divertor_pressure.argsort()
numbers = numbers[arr1inds[::1]]
divertor_pressure = divertor_pressure[arr1inds[::1]]

colormap = cm.cividis
sm = plt.cm.ScalarMappable(cmap=colormap, norm=Normalize(vmin=min(divertor_pressure), vmax=max(divertor_pressure)))

colours = [colormap((P - min(divertor_pressure))/max(divertor_pressure)) for P in divertor_pressure]

filenames_inner = [
    "../data/exposure_conditions_divertor/ITER/{}/{}_inner_target.csv".format(number, number) for number in numbers
    ]
filenames_outer = [
    "../data/exposure_conditions_divertor/ITER/{}/{}_outer_target.csv".format(number, number) for number in numbers
    ]

my_plot = plot_along_divertor(
    filenames=filenames_inner,
    filetypes="ITER",
    quantities=quantities,
    figsize=(6*0.7, 5),
    colors=colours)
for ax in my_plot.axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

my_plot.axs[0].set_ylim(bottom=3e4, top=2e7)
my_plot.axs[1].set_ylim(bottom=1e20, top=3e24)
my_plot.axs[2].set_ylim(bottom=1e-1, top=2e2)


my_plot.axs[-1].get_legend().remove()
my_plot.axs[0].set_ylabel("Heat flux \n (W m$^{-2}$)")
my_plot.axs[1].set_ylabel("Incident flux \n (m$^{-2}$ s$^{-1}$)")
my_plot.axs[2].set_ylabel("Incident energy (eV)")

plt.tight_layout()

my_plot = plot_along_divertor(
    filenames=filenames_outer,
    filetypes="ITER",
    quantities=quantities,
    figsize=(6, 5),
    colors=colours)
for ax in my_plot.axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

my_plot.axs[1].set_ylim(bottom=1e20)
my_plot.axs[2].set_ylim(bottom=1e-1)


custom_lines = [Line2D([0], [0], color="black", linestyle="solid"),
                Line2D([0], [0], color="black", linestyle="dashed")]
my_plot.axs[1].legend(custom_lines, ['Ions', 'Neutrals'], loc="upper right")
my_plot.axs[0].set_ylim(bottom=3e4, top=2e7)
my_plot.axs[1].set_ylim(bottom=1e20, top=3e24)
my_plot.axs[2].set_ylim(bottom=1e-1, top=2e2)
my_plot.axs[-1].get_legend().remove()
my_plot.axs[0].set_ylabel("Heat flux \n (W m$^{-2}$)")
my_plot.axs[1].set_ylabel("Incident flux \n (m$^{-2}$ s$^{-1}$)")
my_plot.axs[2].set_ylabel("Incident energy (eV)")
plt.tight_layout()
plt.colorbar(
    sm, label="Divertor neutral pressure (Pa)",
    ax=my_plot.axs)
plt.show()
