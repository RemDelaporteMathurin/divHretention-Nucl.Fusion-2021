"""
This script has to be executed at the root of the directory
"""
from divHretention import compute_c_max, compute_inventory, \
    plot_along_divertor, Exposition

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np
from scipy.stats import linregress


def as_si(x, ndp):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)

folder = "../data/exposure_conditions_divertor/WEST/"
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

time = 1e7  # s

colormap = cm.inferno
sm = plt.cm.ScalarMappable(cmap=colormap, norm=Normalize(vmin=min(input_powers), vmax=max(input_powers)))
colours = [colormap((IP - min(input_powers))/max(input_powers)) for IP in input_powers]

# inventory at SPs
my_plot = plot_along_divertor(
    filenames=filenames,
    filetypes="WEST",
    quantities=["T_surf", "c_surf", "inventory"],
    figsize=(6, 5),
    colors=colours
)
plt.tight_layout()
plt.colorbar(
    sm, label="Input power (MW)",
    ax=my_plot.axs)

for ax in my_plot.axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
my_plot.axs[-1].get_legend().remove()
plt.savefig("../Figures/WEST/inventory_along_divertor_input_power.pdf")
plt.savefig("../Figures/WEST/inventory_along_divertor_input_power.svg")

# plot integrated inventory in divertor
plt.figure(figsize=(6.4, 2.5))
N_PFU = 480
puff_rates = [2.5e21, 4.44e21]
for puff_rate in puff_rates:
    filenames = [
        folder + "West-LSN-P{:.1e}-IP{:.3f}MW.csv".format(puff_rate, input_power)
        for input_power in input_powers]

    inventories = []
    for filename in filenames:
        res = Exposition(filename, filetype="WEST")
        res.compute_inventory(time)
        inventory = N_PFU*np.trapz(res.inventory, res.arc_length)
        inventories.append(inventory)

    plt.scatter(
        input_powers, inventories, marker="+",
        label="${:s}$".format(as_si(puff_rate, 1)) + " molecule s$^{-1}$")
    res = linregress(np.log10(input_powers), np.log10(inventories))
    a = 10**res.intercept
    b = res.slope
    ip_values = np.logspace(-2, np.log10(input_powers[-1]*1.2))
    line, = plt.plot(
        ip_values, a*ip_values**b,
        linestyle="--")
    plt.annotate(
        "${0:s}".format(as_si(a, 2)) + r"\mathrm{IP} ^{" + "{:.1f}".format(b) + "}$",
        (input_powers[-1]*1.2, a*(input_powers[-1]*1.2)**b),
        color=line.get_color()
    )
plt.xlim(left=0, right=input_powers[-1] + 1)
plt.ylim(bottom=0)
plt.xlabel("Input power (MW)")
plt.ylabel("Divertor H inventory (H)")
plt.legend(loc="lower right")
plt.tight_layout()

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig("../Figures/WEST/inventory_vs_input_power.pdf")
plt.savefig("../Figures/WEST/inventory_vs_input_power.svg")

# plot inventory at SPs
fig, axs = plt.subplots(1, 1, figsize=(6.4, 3), sharey=True, sharex=True)
labels = ["Inner strike point", "Outer strike point"]
linestyles = ["solid", "dashed"]
markers = ["+", "o"]
for i, puff_rate in enumerate([puff_rates[0]]):
    filenames = [
        folder + "West-LSN-P{:.1e}-IP{:.3f}MW.csv".format(puff_rate, input_power)
        for input_power in input_powers]
    inventories_inner_sp, sigmas_inner_sp = [], []
    inventories_outer_sp, sigmas_outer_sp = [], []
    inventories_pz, sigmas_pz = [], []

    for filename in filenames:
        res = Exposition(filename, filetype="WEST")
        res.compute_inventory(time)
        inner_sp_loc_index = np.where(np.abs(res.arc_length-0.20) < 0.005)[0][0]
        outer_sp_loc_index = np.where(np.abs(res.arc_length-0.36) < 0.005)[0][0]
        private_zone_sp_loc_index = np.where(np.abs(res.arc_length-0.28) < 0.005)[0][0]

        inventories_inner_sp.append(res.inventory[inner_sp_loc_index])
        inventories_outer_sp.append(res.inventory[outer_sp_loc_index])
        inventories_pz.append(res.inventory[private_zone_sp_loc_index])

        sigmas_inner_sp.append(res.stdev_inv[inner_sp_loc_index])
        sigmas_outer_sp.append(res.stdev_inv[outer_sp_loc_index])
        sigmas_pz.append(res.stdev_inv[private_zone_sp_loc_index])

    line_spi, = plt.plot(
        input_powers, inventories_inner_sp,
        label="{:.1e}".format(puff_rate) + " molecule s$^{-1}$",
        marker=markers[i],
        color="tab:blue",
        linestyle=linestyles[i])
    line_spo, = plt.plot(
        input_powers, inventories_outer_sp,
        marker=markers[i],
        color="tab:blue",
        linestyle=linestyles[i])
    line_pz,  = plt.plot(
        input_powers, inventories_pz,
        marker=markers[i],
        color="tab:orange",
        linestyle=linestyles[i])

    sigmas_inner_sp = np.array(sigmas_inner_sp)
    sigmas_outer_sp = np.array(sigmas_outer_sp)
    sigmas_pz = np.array(sigmas_pz)

    plt.fill_between(
        input_powers,
        10**(2*sigmas_inner_sp + np.log10(inventories_inner_sp)),
        10**(-2*sigmas_inner_sp + np.log10(inventories_inner_sp)),
        facecolor=line_spi.get_color(), alpha=0.3)
    plt.fill_between(
        input_powers,
        10**(2*sigmas_outer_sp + np.log10(inventories_outer_sp)),
        10**(-2*sigmas_outer_sp + np.log10(inventories_outer_sp)),
        facecolor=line_spo.get_color(), alpha=0.3)
    plt.fill_between(
        input_powers,
        10**(2*sigmas_pz + np.log10(inventories_pz)),
        10**(-2*sigmas_pz + np.log10(inventories_pz)),
        facecolor=line_pz.get_color(), alpha=0.3)

    if i == 0:
        plt.annotate("Inner strike point", (1.05*input_powers[-1], inventories_inner_sp[-1]*1.1), color=line_spi.get_color())
        plt.annotate("Outer strike point", (1.05*input_powers[-1], inventories_outer_sp[-1]*0.9), color=line_spo.get_color())
        plt.annotate("Private zone", (1.05*input_powers[-1], inventories_pz[-1]), color=line_pz.get_color())
plt.ylim(bottom=0, top=1.1e22)
plt.xlim(left=0, right=input_powers[-1] + 0.8)
# plt.legend(loc="upper left")
plt.xlabel("Input power (MW)")
plt.ylabel("Inventory (H m$^{-1}$)")
plt.tight_layout()
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig("../Figures/WEST/inventory_at_sps_and_private_zone_vs_input_power.pdf")
plt.savefig("../Figures/WEST/inventory_at_sps_and_private_zone_vs_input_power.svg")


# #### plot ratios of ions
ratio_ions_inner_sp = []
ratio_ions_outer_sp = []
ratio_ions_private_zone = []

puff_rate = 2.5e21
filenames = [
    folder + "West-LSN-P{:.1e}-IP{:.3f}MW.csv".format(puff_rate, input_power)
    for input_power in input_powers]

for filename in filenames:
    my_exp = Exposition(filename, filetype="WEST")
    res.compute_inventory(time)
    T = 1.1e-4*my_exp.net_heat_flux + 323
    c_max, c_max_ions, c_max_atoms = compute_c_max(
        T,
        my_exp.E_ion, my_exp.E_atom,
        my_exp.angles_ions, my_exp.angles_atoms,
        my_exp.ion_flux, my_exp.atom_flux,
        full_export=True)

    inner_sp_loc_index = np.where(np.abs(res.arc_length-0.20) < 0.005)[0][0]
    outer_sp_loc_index = np.where(np.abs(res.arc_length-0.36) < 0.005)[0][0]
    private_zone_sp_loc_index = np.where(np.abs(res.arc_length-0.28) < 0.005)[0][0]

    ratio_ions_inner_sp.append(c_max_ions[inner_sp_loc_index]/c_max[inner_sp_loc_index])
    ratio_ions_outer_sp.append(c_max_ions[outer_sp_loc_index]/c_max[outer_sp_loc_index])
    ratio_ions_private_zone.append(c_max_ions[private_zone_sp_loc_index]/c_max[private_zone_sp_loc_index])

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(7, 3))
colour_ions = "tab:blue"
colour_atoms = "tab:orange"
line_spo, = axs[0].plot(input_powers, ratio_ions_inner_sp, marker="+", color="tab:blue")
axs[0].fill_between(
    input_powers, np.zeros(len(input_powers)), ratio_ions_inner_sp,
    facecolor=colour_ions, alpha=0.3)
axs[0].fill_between(
    input_powers, np.zeros(len(input_powers)) + 1, ratio_ions_inner_sp,
    facecolor=colour_atoms, alpha=0.3)


axs[0].annotate(r"\textbf{Atoms}", (0.5, 0.8), color=colour_atoms, weight="bold")
axs[0].annotate(r"\textbf{Ions}", (0.6, 0.55), color=colour_ions, weight="bold")

line_spo, = axs[1].plot(input_powers, ratio_ions_outer_sp, marker="+", color="tab:blue")
axs[1].fill_between(
    input_powers, np.zeros(len(input_powers)), ratio_ions_outer_sp,
    facecolor=colour_ions, alpha=0.3)
axs[1].fill_between(
    input_powers, np.zeros(len(input_powers)) + 1, ratio_ions_outer_sp,
    facecolor=colour_atoms, alpha=0.3)


axs[1].annotate(r"\textbf{Atoms}", (0.5, 0.8), color=colour_atoms, weight="bold")
axs[1].annotate(r"\textbf{Ions}", (0.6, 0.5), color=colour_ions, weight="bold")


line_pz, = axs[2].plot(input_powers, ratio_ions_private_zone, marker="+", color="tab:orange")
axs[2].fill_between(
    input_powers, np.zeros(len(input_powers)), ratio_ions_private_zone,
    facecolor=colour_ions, alpha=0.3)
axs[2].fill_between(
    input_powers, np.zeros(len(input_powers)) + 1, ratio_ions_private_zone,
    facecolor=colour_atoms, alpha=0.3)
# axs[2].annotate("Ions", (0.5, 0.1), color="white", weight="bold")
axs[2].annotate(r"\textbf{Atoms}", (0.5, 0.3), color=colour_atoms, weight="bold")

axs[0].set_title("ISP", color=line_spi.get_color())
axs[1].set_title("OSP", color=line_spo.get_color())
axs[2].set_title("Private zone", color=line_pz.get_color())


axs[1].set_xlabel("Input power (MW)")
axs[0].set_ylabel(r"$c_{\mathrm{surface}, \mathrm{ions}} / c_\mathrm{surface}$")

plt.xlim(left=input_powers[0], right=input_powers[-1])
plt.ylim(bottom=0, top=1)
plt.yticks(ticks=[0, 0.5, 1])
plt.tight_layout()
plt.savefig("../Figures/WEST/ions_ratio_vs_input_power.pdf")
plt.savefig("../Figures/WEST/ions_ratio_vs_input_power.svg")

plt.show()
