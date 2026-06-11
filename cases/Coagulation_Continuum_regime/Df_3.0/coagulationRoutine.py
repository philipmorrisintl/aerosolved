#!/usr/bin/python
import numpy as np
import glob
import os
import matplotlib.pyplot as plt

# ===============================
# Input / Output paths
# ===============================
folder = "postProcessing/probes/0/"
out_folder = "combined_output"
os.makedirs(out_folder, exist_ok=True)

# ===============================


# Load simulation parameters
# ===============================
params = np.loadtxt("params.txt")
yMin, yMax = params[8], params[9]       # mass range (kg)
Nsec = int(params[10])                  # number of bins
rhol = params[11]                       # particle density (kg/m^3)
T    = params[5]                        # temperature in K


# ===============================
# Dynamic viscosity using Sutherland's formula
# ===============================
# Constants for air (example values, adjust for your carrier gas)
A_s = 1.67212e-06   # Sutherland constant As
T_s = 170.672       # Sutherland temperature (K)

MU = (A_s * np.sqrt(T)) / (1.0 + (T_s / T))   # dynamic viscosity (kg/m·s)


P=1E5  # pressure
W=28.9596


# ===============================
# Define mass and volume bins
# ===============================
log_mass_edges = np.linspace(np.log10(yMin), np.log10(yMax), Nsec + 1)
mass_edges = 10 ** log_mass_edges
mass_centers = 10 ** ((log_mass_edges[1:] + log_mass_edges[:-1]) / 2)

vol_edges = mass_edges / rhol
vol_centers = mass_centers / rhol
dV = np.diff(vol_edges)

diameters = (6.0 / np.pi * vol_centers) ** (1/3)

# ===============================
# Load rho (time series)
# ===============================
rho_data = np.loadtxt(os.path.join(folder, "rho"))
time = rho_data[:, 0]
rho = rho_data[:, 1]

# ===============================
# Load M and M.* files
# ===============================
M_values = []
headers = ["time", "rho"]

# Load main M file
M_file = os.path.join(folder, "M")
if os.path.isfile(M_file):
    M_main = np.loadtxt(M_file)[:, 1]
    M_values.append(("M", M_main))

# Load additional M.* files
for f in sorted(glob.glob(os.path.join(folder, "M.*"))):
    M = np.loadtxt(f)[:, 1]
    M_values.append((os.path.basename(f), M))

# ===============================
# Compute number concentrations N = M * rho
# ===============================
all_cols = [time, rho]
N_list = []

for name, M in M_values:
    all_cols.append(M)
    headers.append(name)
    
    N = M * rho
    N_name = "N" if name == "M" else name.replace("M", "N", 1)
    all_cols.append(N)
    headers.append(N_name)
    N_list.append(N)

# ===============================
# Compute total particle volume concentration phi
# ===============================
NV_list = []
for idx, (name, _) in enumerate(M_values):
    N = N_list[idx]
    if name == "M":
        NV_list.append(np.zeros_like(N))
    else:
        NV_list.append(N * vol_centers[idx - 1])

phi = np.sum(np.column_stack(NV_list[1:]), axis=1) if len(NV_list) > 1 else np.zeros_like(time)
all_cols.append(phi)
headers.append("phi")

# ===============================
# Compute eta and Si per bin
# ===============================
N_total = N_list[0] if N_list else np.zeros_like(time)
eta_list = []
Si_list = []

for idx, (name, _) in enumerate(M_values):
    if name == "M":
        eta_list.append(np.zeros_like(N_total))
        Si_list.append(np.zeros_like(N_total))
    else:
        bin_idx = idx - 1
        phi_safe = np.where(phi == 0, np.nan, phi)
        eta = (N_total * vol_centers[bin_idx]) / phi_safe
        eta_list.append(eta)

        n_j = np.zeros_like(N_list[idx]) if dV[bin_idx] <= 0 else N_list[idx] / dV[bin_idx]
        N_safe = np.where(N_total == 0, np.nan, N_total)
        Si = (phi * n_j) / (N_safe ** 2)
        Si_list.append(Si)

# Append eta and Si to all_cols for saving
for idx, eta in enumerate(eta_list):
    all_cols.append(eta)
    headers.append(f"eta_bin{idx+1}")
for idx, Si in enumerate(Si_list):
    all_cols.append(Si)
    headers.append(f"Si_bin{idx+1}")

# ===============================
# Save all moments
# ===============================
all_data = np.column_stack(all_cols)
np.savetxt(os.path.join(out_folder, "all_moments.txt"), all_data, fmt="%.6e",
           header="".join(h.ljust(15) for h in headers), comments="")
np.savetxt(os.path.join(out_folder, "all_moments.csv"), all_data, delimiter=",", fmt="%.6e",
           header=",".join(headers), comments="")
print(f"✅ Saved all moments TXT & CSV in {out_folder}")


# ===============================
# Knudsen number & regime
# ===============================
kB = 1.38064852e-23       # J/K
NA = 6.0221409e23         # 1/mol
lambda_mfp = MU / P * np.sqrt(np.pi * kB * T * NA / (2.0 * 0.001 * W))  # m

Kn_list = 2.0 * lambda_mfp / diameters
regime_list = np.array([""]*len(Kn_list), dtype=object)
regime_list[Kn_list > 10] = "Free Molecular"
regime_list[(Kn_list >= 0.1) & (Kn_list <= 10)] = "Transition"
regime_list[Kn_list < 0.1] = "Continuum"

# ===============================
# Plot last-time Si vs eta
# ===============================
last_idx = -1
eta_last = np.array([eta[last_idx] for eta in eta_list])
Si_last  = np.array([Si[last_idx] for Si in Si_list])

plt.figure(figsize=(8,6))
plt.plot(eta_last, Si_last, marker='o', linestyle='-')
plt.xscale('log')
plt.xlabel("eta at final time")
plt.ylabel("Si at final time")
plt.title("Si vs eta for all bins (last time step)")
plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.tight_layout()
plot_file = os.path.join(out_folder, "Si_vs_eta_last_time.png")
plt.savefig(plot_file, dpi=300)
plt.show()
print(f"✅ Saved plot: {plot_file}")

# ===============================
# Save dx and Vx for reference
# ===============================
np.savetxt(os.path.join(out_folder, "dx_Vx.txt"), np.column_stack([diameters, vol_centers]),
           header="dx(m)   Vx(m^3)")

# ===============================
# Save volume-diameter-N-phi-eta-Si-Kn-Regime table (last time step)
# ===============================
outfile_bins = os.path.join(out_folder, "bin_summary.txt")
with open(outfile_bins, "w") as f:
    f.write("{:<15}{:<15}{:<20}{:<15}{:<15}{:<15}{:<15}{:<15}\n".format(
        "volume[m^3]", "diameter[m]", "number_conc[#/m^3]", "phi", "eta", "Si", "Kn", "Regime"))
    for idx, (name, _) in enumerate(M_values):
        if name == "M":  # skip total
            continue
        bin_idx = idx - 1
        f.write("{:<15.6e}{:<15.6e}{:<20.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15.6e}{:<15}\n".format(
            vol_centers[bin_idx],
            diameters[bin_idx],
            N_list[idx][last_idx],
            phi[last_idx],
            eta_list[bin_idx][last_idx],
            Si_list[bin_idx][last_idx],
            Kn_list[bin_idx],
            regime_list[bin_idx]
        ))


# ===============================
# Compute AUC for last-time Si vs eta
# ===============================
auc_last = np.trapz(Si_last, x=eta_last)

# Append AUC to bin_summary.txt
with open(outfile_bins, "a") as f:
    f.write(f"\n# Area under Si vs eta curve at last time = {auc_last:.6e}\n")

print(f"✅ AUC (last time) = {auc_last:.6e} written to {outfile_bins}")


# ===============================
# Plot Si vs eta for multiple times
# ===============================
# List of times you want to plot
target_times = [20000, 20000/2, 20000/5,20000/8]   # <-- edit as needed

plt.figure(figsize=(8,6))

for t in target_times:
    # find closest time index
    time_idx = (np.abs(time - t)).argmin()

    eta_sel = np.array([eta[time_idx] for eta in eta_list])
    Si_sel  = np.array([Si[time_idx] for Si in Si_list])

    plt.plot(eta_sel, Si_sel, marker='o', linestyle='-', label=f"time={time[time_idx]:.3f}")

plt.xscale('log')
plt.xlabel("eta")
plt.ylabel("Si")
plt.title("Si vs eta at selected times")
plt.legend()
plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.tight_layout()

plot_file = os.path.join(out_folder, "Si_vs_eta_multiple_times.png")
plt.savefig(plot_file, dpi=300)
plt.show()
print(f"✅ Saved plot with multiple times: {plot_file}")

# ===============================
# Save Si vs eta side-by-side for multiple times
# ===============================
outfile_side = os.path.join(out_folder, "Si_vs_eta_side_by_side.txt")

# collect data for each target time
data_blocks = []
col_headers = []
for t in target_times:
    # find closest time index
    time_idx = (np.abs(time - t)).argmin()
    eta_sel = np.array([eta[time_idx] for eta in eta_list])
    Si_sel  = np.array([Si[time_idx] for Si in Si_list])

    # stack eta, Si for this time
    block = np.column_stack((eta_sel, Si_sel))
    data_blocks.append(block)

    # headers
    col_headers.append(f"eta_t{int(time[time_idx])}")
    col_headers.append(f"Si_t{int(time[time_idx])}")

# pad shorter arrays if needed (shouldn’t happen normally)
max_len = max(block.shape[0] for block in data_blocks)
for i, block in enumerate(data_blocks):
    if block.shape[0] < max_len:
        pad = np.full((max_len - block.shape[0], 2), np.nan)
        data_blocks[i] = np.vstack([block, pad])

# combine side by side
all_side = np.hstack(data_blocks)

# save to file
header_line = "\t".join(col_headers)
np.savetxt(outfile_side, all_side, fmt="%.6e", header=header_line, comments="")

print(f"✅ Saved side-by-side Si vs eta table: {outfile_side}")




