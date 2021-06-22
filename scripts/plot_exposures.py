import matplotlib.pyplot as plt

from divHretention import *

quantities = [
            "heat_flux",
            "atom_flux", "ion_flux",
            "ion_energy", "atom_energy",
            "ion_angle", "atom_angle",
            ]

filename = "../data/exposure_conditions_divertor/ITER/2398/2398_inner_target.csv"

my_plot = plot_along_divertor(
    filenames=[filename],
    filetypes="ITER",
    quantities=quantities)

for ax in my_plot.axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

filename = "../data/exposure_conditions_divertor/ITER/2398/2398_outer_target.csv"
plot_along_divertor(
    filenames=[filename],
    filetypes="ITER",
    quantities=quantities)

filename = "../data/exposure_conditions_divertor/WEST/West-LSN-P2.5e+21-IP0.449MW.csv"
plot_along_divertor(
    filenames=[filename],
    filetypes="WEST",
    quantities=quantities)

plt.show()
