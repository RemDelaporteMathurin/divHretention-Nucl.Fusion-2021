import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

from scipy.stats import linregress

from divHretention import Exposition


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)

times = np.concatenate((np.logspace(4, 5, num=12, endpoint=False), np.logspace(5, 7, num=10)))
time_in_shots = times/400

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(6.4, 5))
folder = "../data/exposure_conditions_divertor/ITER/"

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


# sort arrays
divertor_pressure = np.array(divertor_pressure)
numbers = np.array(numbers)
arr1inds = divertor_pressure.argsort()
numbers = numbers[arr1inds[::1]]
divertor_pressure = divertor_pressure[arr1inds[::1]]

colormap = cm.cividis
sm = plt.cm.ScalarMappable(cmap=colormap, norm=Normalize(vmin=min(divertor_pressure), vmax=max(divertor_pressure)))

filenames_inner = [
    "../data/exposure_conditions_divertor/ITER/{}/{}_inner_target.csv".format(number, number) for number in numbers
    ]
filenames_outer = [
    "../data/exposure_conditions_divertor/ITER/{}/{}_outer_target.csv".format(number, number) for number in numbers
    ]
N_cassettes = 54
N_IVT = N_cassettes*16
N_OVT = N_cassettes*22
avogadro = 6.022e23
as_, bs_ = [], []
for file_inner, file_outer, pressure in zip(filenames_inner, filenames_outer, divertor_pressure):
    inventories = []
    for time in times:
        inventory = 0
        for filename, N in zip([file_inner, file_outer], [N_IVT, N_OVT]):
            res = Exposition(filename, filetype="ITER")
            res.compute_inventory(time)
            inventory += N*np.trapz(res.inventory, res.arc_length)
        inventories.append(inventory)
    inventories = np.array(inventories)/avogadro  # in g
    plt.sca(axs[0])
    plt.plot(time_in_shots, inventories, marker="+", color=colormap((pressure - min(divertor_pressure))/max(divertor_pressure)))
    plt.sca(axs[1])
    inv_dt = np.diff(inventories)/np.diff(time_in_shots)
    inv_dt *= 1e3  # in mg
    plt.plot(
        time_in_shots[:-1], inv_dt,
        color=colormap((pressure - min(divertor_pressure))/max(divertor_pressure)))
plt.sca(axs[0])
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Divertor \n inventory (g)")


# plot additional inventory
plt.sca(axs[1])

plt.ylim(bottom=0)
plt.xscale("log")
plt.xlabel("Number of discharges")
plt.ylabel("Additionnal inventory \n per shot (mg/shot)")

plt.colorbar(sm, label="Divertor neutral pressure (Pa)", ax=axs)
# plt.tight_layout()

for ax in axs:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
plt.savefig("../Figures/ITER/inventory_vs_time.pdf")
plt.savefig("../Figures/ITER/inventory_vs_time.svg")
plt.show()
