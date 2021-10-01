"""
Figures 9, 10, 11
"""
from divHretention import compute_c_max, compute_inventory, \
    plot_along_divertor, Exposition
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.stats import linregress

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)


def as_si(x, ndp):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))


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

N_PFU = 480
time = 1e7  # s
# #### plot inventory along divertor

my_plot = plot_along_divertor(
    filenames,
    filetypes="WEST",
    quantities=["T_surf", "c_surf", "inventory"],
    colors=colours,
    figsize=(6, 5))

plt.tight_layout()
plt.colorbar(
    sm, label="Puffing rate (molecule s$^{-1}$)",
    ax=my_plot.axs)

for ax in my_plot.axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

my_plot.axs[-1].get_legend().remove()
plt.savefig('../Figures/WEST/inventory_along_divertor.pdf')
plt.savefig('../Figures/WEST/inventory_along_divertor.svg')

# #### plot contribution ions along divertor
inventory_strike_point_inner = []
inventory_strike_point_outer = []
inventory_private_zone = []

sigma_strike_point_inner = []
sigma_strike_point_outer = []
sigma_private_zone = []

ratio_ions_inner_sp = []
ratio_ions_outer_sp = []
ratio_ions_private_zone = []

integrated_inventories = []
plt.figure()
for i, filename in enumerate(filenames):
    my_exp = Exposition(filename, filetype="WEST")
    arc_length = my_exp.arc_length
    T = 1.1e-4*my_exp.net_heat_flux + 323
    c_max, c_max_ions, c_max_atoms = compute_c_max(
        T,
        my_exp.E_ion, my_exp.E_atom,
        my_exp.angles_ions, my_exp.angles_atoms,
        my_exp.ion_flux, my_exp.atom_flux,
        full_export=True)
    inventories, sigmas = compute_inventory(T, c_max, time=time)
    integrated_inventories.append(N_PFU*np.trapz(inventories, arc_length))

    line, = plt.plot(arc_length, c_max_ions/c_max, color=colours[i])

    inner_sp_loc_index = np.where(np.abs(arc_length-0.20) < 0.005)[0][0]
    outer_sp_loc_index = np.where(np.abs(arc_length-0.36) < 0.005)[0][0]
    private_zone_sp_loc_index = np.where(np.abs(arc_length-0.28) < 0.005)[0][0]

    inventory_strike_point_inner.append(inventories[inner_sp_loc_index])
    inventory_strike_point_outer.append(inventories[outer_sp_loc_index])
    inventory_private_zone.append(inventories[private_zone_sp_loc_index])

    sigma_strike_point_inner.append(sigmas[inner_sp_loc_index])
    sigma_strike_point_outer.append(sigmas[outer_sp_loc_index])
    sigma_private_zone.append(sigmas[private_zone_sp_loc_index])

    ratio_ions_inner_sp.append(c_max_ions[inner_sp_loc_index]/c_max[inner_sp_loc_index])
    ratio_ions_outer_sp.append(c_max_ions[outer_sp_loc_index]/c_max[outer_sp_loc_index])
    ratio_ions_private_zone.append(c_max_ions[private_zone_sp_loc_index]/c_max[private_zone_sp_loc_index])

plt.colorbar(sm, label="Puffing rate (molecule s$^{-1}$)")
plt.ylim(0, 1)
plt.xlabel("Distance along divertor (m)")
plt.ylabel("c surface (ions) / c surface")
plt.savefig('../Figures/WEST/ion_ratio_along_divertor.pdf')
plt.savefig('../Figures/WEST/ion_ratio_along_divertor.svg')

# #### plot inventory vs puffing rate

res = linregress(np.log10(Ps), np.log10(integrated_inventories))

puffin_rate_values = np.linspace(0, 5e21, num=200)

plt.figure(figsize=(6.4, 2.5))
line, = plt.plot(puffin_rate_values, 10**res.intercept*puffin_rate_values**res.slope, linestyle="--")
plt.annotate(
    "${0:s}".format(as_si(10**res.intercept, 2)) + r"p ^ {" + "{:.1f}".format(res.slope) + "}$",
    (3e21, 4.4e23),
    color=line.get_color())
plt.scatter(Ps, integrated_inventories, marker="+")
plt.xlabel("Puffing rate (molecule s$^{-1}$)")
plt.ylabel("Divertor H inventory (H)")
plt.ylim(bottom=0)
plt.xlim(left=0)

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('../Figures/WEST/inventory_vs_puffing_rate.pdf')
plt.savefig('../Figures/WEST/inventory_vs_puffing_rate.svg')

# #### plot inventory at SPs and private zone
plt.figure(figsize=(6.4, 3))
line_spi, = plt.plot(Ps, inventory_strike_point_inner, marker="+", color="tab:blue")
line_spo, = plt.plot(Ps, inventory_strike_point_outer, marker="+", color="tab:blue")
line_pz,  = plt.plot(Ps, inventory_private_zone, marker="+", color="tab:orange")

sigma_strike_point_inner = np.array(sigma_strike_point_inner)
sigma_strike_point_outer = np.array(sigma_strike_point_outer)
sigma_private_zone = np.array(sigma_private_zone)

plt.fill_between(
    Ps,
    10**(2*sigma_strike_point_inner + np.log10(inventory_strike_point_inner)),
    10**(-2*sigma_strike_point_inner + np.log10(inventory_strike_point_inner)),
    facecolor=line_spi.get_color(), alpha=0.3)
plt.fill_between(
    Ps,
    10**(2*sigma_strike_point_outer + np.log10(inventory_strike_point_outer)),
    10**(-2*sigma_strike_point_outer + np.log10(inventory_strike_point_outer)),
    facecolor=line_spo.get_color(), alpha=0.3)
plt.fill_between(
    Ps,
    10**(2*sigma_private_zone + np.log10(inventory_private_zone)),
    10**(-2*sigma_private_zone + np.log10(inventory_private_zone)),
    facecolor=line_pz.get_color(), alpha=0.3)

plt.annotate(
    "Inner strike point",
    (1.05*Ps[-1], 0.91*inventory_strike_point_inner[-1]),
    color=line_spi.get_color())
plt.annotate(
    "Outer strike point",
    (1.05*Ps[-1], inventory_strike_point_outer[-1]),
    color=line_spo.get_color())
plt.annotate("Private zone", (1.05*Ps[-1], inventory_private_zone[-1]), color=line_pz.get_color())
plt.xlim(right=1.4*Ps[-1])
plt.xlabel("Puffing rate (molecule s$^{-1}$)")
plt.ylabel("Inventory (H/m)")
plt.ylim(bottom=0)
plt.xlim(left=0)
plt.tight_layout()
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig('../Figures/WEST/inventory_at_sp_and_private_zone.pdf')
plt.savefig('../Figures/WEST/inventory_at_sp_and_private_zone.svg')

# #### plot ions vs atoms
fig, axs = plt.subplots(1, 3, sharey="row", sharex=True, figsize=(7, 3))
line_spi, = axs[0].plot(Ps, ratio_ions_inner_sp, marker="+", color="tab:blue")
axs[0].fill_between(
    Ps, np.zeros(len(Ps)), ratio_ions_inner_sp,
    facecolor='tab:blue', alpha=0.3)
axs[0].fill_between(
    Ps, np.zeros(len(Ps)) + 1, ratio_ions_inner_sp,
    facecolor='tab:orange', alpha=0.3)

line_spo, = axs[1].plot(Ps, ratio_ions_outer_sp, marker="+", color="tab:blue")
axs[1].fill_between(
    Ps, np.zeros(len(Ps)), ratio_ions_outer_sp,
    facecolor='tab:blue', alpha=0.3)
axs[1].fill_between(
    Ps, np.zeros(len(Ps)) + 1, ratio_ions_outer_sp,
    facecolor='tab:orange', alpha=0.3)

line_pz,  = axs[2].plot(Ps, ratio_ions_private_zone, marker="+", color="tab:orange")
axs[2].fill_between(
    Ps, np.zeros(len(Ps)), ratio_ions_private_zone,
    facecolor='tab:blue', alpha=0.3)
axs[2].fill_between(
    Ps, np.zeros(len(Ps)) + 1, ratio_ions_private_zone,
    facecolor='tab:orange', alpha=0.3)

axs[0].set_title("ISP", color=line_spi.get_color())
axs[1].set_title("OSP", color=line_spo.get_color())
axs[2].set_title("Private zone", color=line_pz.get_color())

colour_ions = "tab:blue"
colour_atoms = "tab:orange"
axs[0].annotate(r"\textbf{Ions}", (3e21, 0.5), color=colour_ions, weight="bold")
axs[0].annotate(r"\textbf{Atoms}", (3.4e21, 0.6), color=colour_atoms, weight="bold")
axs[1].annotate(r"\textbf{Ions}", (3e21, 0.4), color=colour_ions, weight="bold")
axs[1].annotate(r"\textbf{Atoms}", (3.4e21, 0.5), color=colour_atoms, weight="bold")

# axs[2].annotate("Ions", (0.5e21, 0.1), color="grey", weight="bold")
axs[2].annotate(r"\textbf{Atoms}", (1e21, 0.5), color=colour_atoms, weight="bold")

plt.sca(axs[0])
plt.xlim(left=Ps[0], right=Ps[-1])
axs[1].set_xlabel("Puffing rate \n (molecule s$^{-1}$)")
plt.ylabel(r"$c_{\mathrm{surface}, \mathrm{ions}} / c_\mathrm{surface}$")
plt.ylim(bottom=0, top=1)
plt.yticks(ticks=[0, 0.5, 1])
# plt.xticks(ticks=[Ps[0], 2e21, Ps[-1]])
plt.tight_layout()
plt.savefig('../Figures/WEST/ion_ratio_at_sp_and_private_zone.pdf')
plt.savefig('../Figures/WEST/ion_ratio_at_sp_and_private_zone.svg')

plt.show()
