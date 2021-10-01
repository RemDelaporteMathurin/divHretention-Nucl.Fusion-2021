"""
Figures 13, 14, 15
"""

import re
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

import numpy as np
from divHretention import plot_along_divertor, \
    Exposition, compute_c_max
import divHretention
from scipy.interpolate import interp1d

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)

divHretention.step_mb = 1

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

time = 1e7  # s
# sort arrays
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

res_inner, res_outer = [], []

for filename in filenames_inner:
    res = Exposition(filename, filetype="ITER")
    res.compute_inventory(time)
    res_inner.append(res)

for i, filename in enumerate(filenames_outer):
    res = Exposition(filename, filetype="ITER")
    res.compute_inventory(time)
    res_outer.append(res)

# plot particle fluxes
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=[6.4, 4])

for i, res in enumerate(res_inner):
    flux_atm = res.atom_flux
    flux_ion = res.ion_flux
    arc_length = res.arc_length
    ax1.plot(arc_length, flux_ion, color=colours[i])
    ax2.plot(arc_length, flux_atm, color=colours[i])

ax1.annotate("Ions", (0.6, 1e24))
ax2.annotate("Atoms", (0.6, 1e24))

ax1.set_xlabel("Distance along divertor (m)")
ax2.set_xlabel("Distance along divertor (m)")
ax1.set_ylabel("Particle flux (m$^{-2}$.s$^{-1}$)")
plt.yscale("log")
plt.tight_layout()
plt.colorbar(
    sm, label="Divertor neutral pressure (Pa)",
    ax=ax2)
plt.savefig('../Figures/ITER/fluxes_distribution.svg')
plt.savefig('../Figures/ITER/fluxes_distribution.pdf')

# plot integrated inventory
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=[6.4, 4])

ax1.set_title("IVT")
ax2.set_title("OVT")

IVT_integrals_ion, IVT_integrals_atom = [], []
for i, res in enumerate(res_inner):
    flux_atm = res.atom_flux
    flux_ion = res.ion_flux
    arc_length = res.arc_length
    IVT_integrals_ion.append(np.trapz(flux_ion, arc_length))
    IVT_integrals_atom.append(np.trapz(flux_atm, arc_length))

ax1.plot(divertor_pressure, IVT_integrals_ion)
ax1.plot(divertor_pressure, IVT_integrals_atom, linestyle="dashed")

OVT_integrals_ion, OVT_integrals_atom = [], []

for i, res in enumerate(res_outer):
    flux_atm = res.atom_flux
    flux_ion = res.ion_flux
    arc_length = res.arc_length
    OVT_integrals_ion.append(np.trapz(flux_ion, arc_length))
    OVT_integrals_atom.append(np.trapz(flux_atm, arc_length))
ax2.plot(divertor_pressure, OVT_integrals_ion)
ax2.plot(divertor_pressure, OVT_integrals_atom, linestyle="dashed")

ax1.set_ylabel(r"$\int \varphi$ (m$^{-1}$ s$^{-1}$)")
ax1.set_xlabel("Divertor neutral pressure (Pa)")
ax2.set_xlabel("Divertor neutral pressure (Pa)")

ax2.annotate(
    "Ions", (divertor_pressure[-1], OVT_integrals_ion[-1]),
    color="tab:blue")
ax2.annotate(
    "Atoms", (divertor_pressure[-1], OVT_integrals_atom[-1]),
    color="tab:orange")
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
plt.ylim(bottom=0)
plt.tight_layout()
plt.savefig('../Figures/ITER/integral_fluxes.svg')
plt.savefig('../Figures/ITER/integral_fluxes.pdf')

# compute total inventory
def integrate_inventory(res, x_max='max'):
    inventory_interp = interp1d(res.arc_length, res.inventory)
    if x_max == 'max':
        x_limit = res.arc_length.max()
    else:
        x_limit = x_max
    x_values = np.linspace(res.arc_length.min(), x_limit)
    inventory = np.trapz(inventory_interp(x_values), x_values)
    return inventory


def compute_total_inventory(x_max='max'):
    inventories_IVT, inventories_OVT = [], []
    for res in res_inner:
        inventory = integrate_inventory(res)
        inventories_IVT.append(inventory)
    for res in res_outer:
        inventory = integrate_inventory(res)
        inventories_OVT.append(inventory)
    return inventories_IVT, inventories_OVT


# ###### plot inventory vs divertor pressure
N_cassettes = 54
N_IVT = N_cassettes*16
N_OVT = N_cassettes*22

inventories_IVT, inventories_OVT = compute_total_inventory()

avogadro = 6.022e23
inventories_IVT, inventories_OVT = \
    N_IVT/avogadro*np.array(inventories_IVT), \
    N_OVT/avogadro*np.array(inventories_OVT)
inventories_total = inventories_IVT + inventories_OVT
plt.figure()
line_tot, = plt.plot(divertor_pressure, inventories_total, marker="+", color="tab:blue")
line_inner, = plt.plot(divertor_pressure, inventories_IVT, marker="+", color="tab:blue")
plt.fill_between(
    divertor_pressure, np.zeros(len(divertor_pressure)), inventories_IVT,
    alpha=0.6, color=line_inner.get_color())
plt.fill_between(
    divertor_pressure, inventories_IVT, inventories_total,
    alpha=0.3, color=line_tot.get_color())
plt.xlabel("Divertor neutral pressure (Pa)")
plt.ylabel("Divertor H inventory (g)")

plt.ylim(bottom=0)
plt.xlim(left=divertor_pressure[0], right=divertor_pressure[-1] + 1.5)
plt.annotate(r"\textbf{Inner Target}", (7, 2), color="tab:blue", weight="bold")
plt.annotate(r"\textbf{Outer Target}", (7, 10), color="tab:blue", weight="bold")

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig('../Figures/ITER/inventory_vs_divertor_pressure.pdf')
plt.savefig('../Figures/ITER/inventory_vs_divertor_pressure.svg')

# ###### plot exposure conditions along divertor

my_plot_inner = plot_along_divertor(
    filenames=filenames_inner,
    filetypes="ITER",
    quantities=["T_surf", "c_surf", "inventory"],
    colors=colours,
    figsize=(6*0.7, 5))
plt.tight_layout()
# my_plot_inner.axs[0].annotate("IVT", (0.5, 800))
# plt.colorbar(
#     sm, label="Divertor neutral pressure (Pa)",
#     ax=my_plot_inner.axs)

for ax in my_plot_inner.axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

my_plot_inner.axs[0].set_ylim(top=2100)
my_plot_inner.axs[1].set_ylim(6e19, 7e22)
my_plot_inner.axs[2].set_ylim(1e20, 1e22)
my_plot_inner.axs[-1].get_legend().remove()
plt.savefig('../Figures/ITER/inventory_along_inner_divertor.pdf')
plt.savefig('../Figures/ITER/inventory_along_inner_divertor.svg')

my_plot_outer = plot_along_divertor(
    filenames=filenames_outer,
    filetypes="ITER",
    quantities=["T_surf", "c_surf", "inventory"],
    colors=colours,
    figsize=(6, 5))
plt.tight_layout()
# my_plot_outer.axs[0].annotate("OVT", (0.5, 1500))
plt.colorbar(
    sm, label="Divertor neutral pressure (Pa)",
    ax=my_plot_outer.axs)

for ax in my_plot_outer.axs:
    # ax.set_ylabel("")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

my_plot_outer.axs[0].set_ylim(top=2100)
my_plot_outer.axs[1].set_ylim(6e19, 7e22)
my_plot_outer.axs[2].set_ylim(1e20, 1e22)

my_plot_outer.axs[-1].get_legend().remove()
plt.savefig('../Figures/ITER/inventory_along_outer_divertor.pdf')
plt.savefig('../Figures/ITER/inventory_along_outer_divertor.svg')

# ###### compute inventory as function of x
fig, axs = plt.subplots(1, 2, sharey=True, figsize=(6.4, 3))
titles = ["IVT", "OVT"]
for i, results in enumerate([res_inner, res_outer]):
    plt.sca(axs[i])
    plt.title(titles[i])
    plt.xlabel("Distance along divertor (m)")
    for j, res in enumerate(results):
        inventories = []
        x_max_values = np.linspace(res.arc_length.min(), res.arc_length.max())
        for x_max in x_max_values:
            inventories.append(integrate_inventory(res, x_max=x_max))
        inventories = np.array(inventories)
        inventories *= 1/inventories.max()
        plt.plot(x_max_values, inventories, color=colours[j], alpha=0.8)

axs[0].set_ylabel("Cumulative inventory (norm.)")
plt.yticks(ticks=[0, 0.5, 1])
fig.colorbar(sm, label="Divertor neutral pressure (Pa)")
plt.tight_layout()
plt.subplots_adjust(wspace=0)
plt.ylim(0, 1)

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig('../Figures/ITER/cumulative_inventory.pdf')
plt.savefig('../Figures/ITER/cumulative_inventory.svg')

# ###### plot inventory at SPs
plt.figure(figsize=(6.4, 3))
labels = ["Inner strike point", "Outer strike point"]
for i, results in enumerate([res_inner, res_outer]):
    inventories, sigmas = [], []
    for res in results:
        sp_loc_index = np.where(res.temperature == res.temperature.max())[0][0]
        inventories.append(res.inventory[sp_loc_index])
        sigmas.append(res.stdev_inv[sp_loc_index])
    sigmas = np.array(sigmas)
    line, = plt.plot(divertor_pressure, inventories, label=labels[i], marker="+", color="tab:blue")
    plt.fill_between(
        divertor_pressure,
        10**(2*sigmas + np.log10(inventories)),
        10**(-2*sigmas + np.log10(inventories)),
        facecolor=line.get_color(), alpha=0.3)
    plt.annotate(labels[i], (1.05*divertor_pressure[-1], inventories[-1]), color=line.get_color())

plt.xlabel("Divertor neutral pressure (Pa)")
plt.ylabel("H inventory (H m$^{-1}$)")
plt.xlim(left=0, right=divertor_pressure[-1] + 4.5)
plt.ylim(bottom=0)
plt.tight_layout()

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig('../Figures/ITER/inventory_at_strike_points.pdf')
plt.savefig('../Figures/ITER/inventory_at_strike_points.svg')

# ###### plot neutral contribution
ratios = [[], []]
for i, results in enumerate([filenames_inner, filenames_outer]):
    for filename in results:
        my_exp = Exposition(filename, filetype="ITER")

        T = 1.1e-4*my_exp.net_heat_flux + 323
        c_max, c_max_ions, c_max_atoms = compute_c_max(
            T,
            my_exp.E_ion, my_exp.E_atom,
            my_exp.angles_ions, my_exp.angles_atoms,
            my_exp.ion_flux, my_exp.atom_flux,
            full_export=True)

        inner_sp_loc_index = np.where(res.temperature == res.temperature.max())[0][0]

        ratios[i].append(c_max_ions[inner_sp_loc_index]/c_max[inner_sp_loc_index])

fig, axs = plt.subplots(1, 2, sharey="row", sharex=True, figsize=(5.5, 3))
colour_ions = "tab:blue"
colour_atoms = "tab:orange"
ratio_ions_inner_sp = ratios[0]
line_spi, = axs[0].plot(divertor_pressure, ratio_ions_inner_sp, marker="+", color="tab:blue")
axs[0].fill_between(
    divertor_pressure, np.zeros(len(divertor_pressure)), ratio_ions_inner_sp,
    facecolor=colour_ions, alpha=0.3)
axs[0].fill_between(
    divertor_pressure, np.zeros(len(divertor_pressure)) + 1, ratio_ions_inner_sp,
    facecolor=colour_atoms, alpha=0.3)

ratio_ions_outer_sp = ratios[1]
line_spo, = axs[1].plot(divertor_pressure, ratio_ions_outer_sp, marker="+", color="tab:blue")
axs[1].fill_between(
    divertor_pressure, np.zeros(len(divertor_pressure)), ratio_ions_outer_sp,
    facecolor=colour_ions, alpha=0.3)
axs[1].fill_between(
    divertor_pressure, np.zeros(len(divertor_pressure)) + 1, ratio_ions_outer_sp,
    facecolor=colour_atoms, alpha=0.3)

axs[0].set_title("ISP")
axs[1].set_title("OSP")

axs[0].annotate(r"\textbf{Ions}", (3, 0.5), color=colour_ions, weight="bold")
axs[0].annotate(r"\textbf{Atoms}", (3.4, 0.7), color=colour_atoms, weight="bold")
axs[1].annotate(r"\textbf{Ions}", (3.4, 0.55), color=colour_ions, weight="bold")
axs[1].annotate(r"\textbf{Atoms}", (5, 0.7), color=colour_atoms, weight="bold")


plt.sca(axs[0])
plt.xlim(left=divertor_pressure[0], right=divertor_pressure[-1])
axs[0].set_xlabel("Divertor neutral pressure (Pa)")
axs[1].set_xlabel("Divertor neutral pressure (Pa)")
plt.ylabel(r"$c_{\mathrm{surface}, \mathrm{ions}} / c_\mathrm{surface}$")
plt.ylim(bottom=0, top=1)
plt.yticks(ticks=[0, 0.5, 1])
plt.tight_layout()
plt.savefig('../Figures/ITER/ratio_ions_atoms.pdf')
plt.savefig('../Figures/ITER/ratio_ions_atoms.svg')

plt.show()
