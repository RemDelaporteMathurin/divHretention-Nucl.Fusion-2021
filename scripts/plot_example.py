import matplotlib.pyplot as plt
from divHretention import plot_Tc_map_with_subplots, DEFAULT_TIME, \
    compute_inventory


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)

filename = "../data/exposure_conditions_divertor/WEST/West-LSN-P2.5e+21-IP0.449MW.csv"


# trigger DEFAULT_TIME
# compute_inventory([300], [1e20], DEFAULT_TIME)

my_plot = plot_Tc_map_with_subplots(
    filenames=[filename], filetypes="WEST",
    T_bounds=[320, 400], figsize=(5, 5))
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
my_plot.show()
