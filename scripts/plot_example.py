"""
Figure 3
"""
import matplotlib.pyplot as plt
from divHretention import plot_Tc_map_with_subplots, compute_inventory
from matplotlib import cm
from matplotlib.colors import LogNorm


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)

# ITER
filename = "../data/exposure_conditions_divertor/ITER/2398/2398_inner_target.csv"

my_plot = plot_Tc_map_with_subplots(
    filenames=[filename], filetypes="ITER",
    T_bounds=[320, 1200], figsize=(5, 5))
my_plot.axs_bottom[1].get_yaxis().set_ticks([])
my_plot.axs_top[1].yaxis.set_label_position("right")
my_plot.axs_top[1].yaxis.tick_right()
my_plot.axs_top[1].set_ylabel("Inventory per unit \n thickness (H/m)")

# change everything to orange
for ax in [*my_plot.axs_top, *my_plot.axs_bottom]:
    for line in ax.get_lines():
        line.set_color("tab:orange")

for cl in my_plot.axs_top[1].collections:
    cl.set_facecolor("tab:orange")

plt.tight_layout()


# WEST
filename = "../data/exposure_conditions_divertor/WEST/West-LSN-P2.5e+21-IP0.449MW.csv"

my_plot = plot_Tc_map_with_subplots(
    filenames=[filename], filetypes="WEST",
    T_bounds=[320, 400], figsize=(5, 5))
my_plot.axs_top[1].set_ylabel("Inventory (H/m)")
for ax in [*my_plot.axs_top, my_plot.axs_bottom[1]]:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# change everything to orange
for ax in [*my_plot.axs_top, *my_plot.axs_bottom]:
    for line in ax.get_lines():
        line.set_color("tab:orange")

for cl in my_plot.axs_top[1].collections:
    cl.set_facecolor("tab:orange")

plt.tight_layout()

colormap = cm.viridis
inv_min = compute_inventory([320], [1e20], 1e7)[0][0]
inv_max = compute_inventory([400], [1e23], 1e7)[0][0]
sm = plt.cm.ScalarMappable(cmap=colormap, norm=LogNorm(vmin=inv_min, vmax=inv_max))

plt.colorbar(
    sm, label="Inventory per unit thickness (H/m)",
    ax=[*my_plot.axs_bottom, *my_plot.axs_top])

my_plot.show()
