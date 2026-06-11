#!/usr/bin/python

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.special as sps

# Add custom scripts folder to path
sys.path.insert(0, '../../scripts')
#import figStyle as fs

# ===============================
# Settings
# ===============================
folders = [
    ("Df_3.0/1E9",  "Df =3.0"),
    ("Df_1.8/1E9",  "Df =1.8"),
]

colors = ['#1f77b4', '#ff7f0e', 'green', 'purple', 'brown', 'pink']

# Output folder
out_folder = "output"
os.makedirs(out_folder, exist_ok=True)

# ===============================
# Functions
# ===============================
def load_params(folder):
    data = np.loadtxt(os.path.join(folder, "params.txt"))
    return {
        "Kn0": data[0],
        "K0": data[1],
        "Kfm": data[2],
        "M0": data[3],
        "tau": data[4],
        "T": data[5],
        "l": data[6],
        "sigma0": data[7],
        "yMin": data[8],
        "yMax": data[9],
        "Nsec": int(data[10]),
        "rhol": data[11]
    }

def load_rho_time(folder):
    rho_time = np.loadtxt(os.path.join(folder, "postProcessing/probes/0/rho"))
    return rho_time[:,0], rho_time[:,1]

def compute_bins(yMin, yMax, Nsec, rhol):
    logy = np.linspace(np.log10(yMin), np.log10(yMax), Nsec+1)
    y = 10**logy
    x = 10**((logy[1:] + logy[:-1]) / 2.0)
    dx = (x / rhol * 6.0 / np.pi)**(1/3) * 1e6
    dy = (y / rhol * 6.0 / np.pi)**(1/3) * 1e6
    return x, dx, y, dy

def compute_dn_dlogd(folder, Nsec, rho, idx, dy, rhol):
    N = np.zeros(Nsec)
    for j in range(1, Nsec):
        fileName = f'M.{str(j).zfill(int(np.log10(Nsec))+1)}'
        M = np.loadtxt(os.path.join(folder, "postProcessing/probes/0", fileName))[idx, 1]
        N[j] = M * rho[idx]
    dlogd_sec = np.log10(dy[1:Nsec+1]) - np.log10(dy[0:Nsec])
    sim_dnldlogd = N / dlogd_sec / 1e6
    return N, sim_dnldlogd

def save_txt(file_path, data, header):
    np.savetxt(file_path, data, fmt="%.6e", header=header, comments="")

# ===============================
# Prepare Figures
# ===============================
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

combined_dn_dlogdp = {}
combined_N_by_N0 = {}

# Initial distribution point
data0 = np.loadtxt(os.path.join(folders[0][0], "params.txt"))
M0_init = data0[3]
rhog = 1.2235553
d0_um = 1000e-9 * 1e6
N0_cm3 = M0_init * rhog / 1e6
ax1.plot(d0_um, N0_cm3, 'o', color='black', label="Initial (t = 0 h)")
ax1.axvline(d0_um, linestyle=":", color='black', linewidth=1)

# ===============================
# Process all folders
# ===============================
for (folder, label), color in zip(folders, colors):
    params = load_params(folder)
    time, rho = load_rho_time(folder)
    x, dx, y, dy = compute_bins(params["yMin"], params["yMax"], params["Nsec"], params["rhol"])
    
    idx = -1
    N, sim_dnldlogd = compute_dn_dlogd(folder, params["Nsec"], rho, idx, dy, params["rhol"])
    
    mask_sim = sim_dnldlogd > 0
    ax1.plot(dx[mask_sim], sim_dnldlogd[mask_sim], label=label, color=color)
    
    sim_data = np.column_stack((dx, N / 1e6, sim_dnldlogd))
    save_txt(os.path.join(out_folder, f"{label.replace(' ', '_')}_dn_dlogdp.txt"),
             sim_data, "d[um]   n[#/cm3]   dn/dlog10(Dp)[#/cm3]")

    Mfile = np.loadtxt(os.path.join(folder, "postProcessing/probes/0/M"))
    t_sim = Mfile[:,0] / params["tau"]
    N_by_N0 = Mfile[:,1] / params["M0"]
    ax2.plot(t_sim, N_by_N0, label=label, color=color)
    
    nn0_data = np.column_stack((t_sim, N_by_N0))
    save_txt(os.path.join(out_folder, f"{label.replace(' ', '_')}_N_by_N0.txt"),
             nn0_data, "t/tau   N/N0")
    
    combined_dn_dlogdp[label] = np.column_stack((dx[mask_sim], sim_dnldlogd[mask_sim]))
    combined_N_by_N0[label] = np.column_stack((t_sim, N_by_N0))

# ===============================
# Park simplified
# ===============================
params0 = load_params(folders[0][0])

Z0 = np.log(params0["sigma0"])**2
p = 1.0 + np.exp(Z0)

t_all = np.concatenate([combined_N_by_N0[label][:,0] for label in combined_N_by_N0])
t_min, t_max = np.min(t_all), np.max(t_all)
t_common = np.logspace(np.log10(t_min), np.log10(t_max), 300)

N_park = 1.0 / (1.0 + p * t_common)

ax2.plot(t_common, N_park, '--k', linewidth=2, label='Park (simplified)')

# Store Park
combined_N_by_N0["Park (simplified)"] = np.column_stack((t_common, N_park))

# ===============================
# Save combined data
# ===============================
def save_combined_N_by_N0(path, data_dict, folders):
    all_labels = [label for _, label in folders] + ["Park (simplified)"]

    tmin = min(np.min(arr[:,0]) for arr in data_dict.values())
    tmax = max(np.max(arr[:,0]) for arr in data_dict.values())
    common_t = np.logspace(np.log10(tmin), np.log10(tmax), 200)

    with open(path, "w") as f:
        f.write("t/tau " + " ".join(all_labels) + "\n")

        for t in common_t:
            row = [f"{t:.6e}"]
            for label in all_labels:
                arr = data_dict[label]
                row.append(f"{np.interp(t, arr[:,0], arr[:,1]):.6e}")
            f.write(" ".join(row) + "\n")

save_combined_N_by_N0(os.path.join(out_folder, "combined_N_by_N0.txt"),
                     combined_N_by_N0, folders)

# ===============================
# Style and Save Plots
# ===============================
def style_plot(ax, xlabel, ylabel):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(xlabel, fontsize=14, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=7, frameon=False)

style_plot(ax1, r'Particle diameter $D_p$ [$\mu$m]', r'$dn/d\log_{10}D_p$')
fig1.tight_layout()
fig1.savefig(os.path.join(out_folder, 'comparison_dn_dlogdp.png'), dpi=300)

style_plot(ax2, r'$t/\tau$', r'$N/N_0$')
fig2.tight_layout()
fig2.savefig(os.path.join(out_folder, 'comparison_N_by_N0.png'), dpi=300)
