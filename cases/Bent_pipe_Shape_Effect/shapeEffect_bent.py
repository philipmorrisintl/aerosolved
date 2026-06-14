import numpy as np
import matplotlib.pyplot as plt

# Define the folder paths with legend-friendly labels
folders = {
    'Shape factor = 1': 'Shape_1/postProcessing/dropletFlux/0',
    'Shape factor = 1.5': 'Shape_1.5/postProcessing/dropletFlux/0'
}

# Define particle sizes (same for both cases)
particle_sizes = np.array([0.06, 0.07, 0.09, 0.12, 0.15, 0.19, 0.25, 0.32, 0.41, 0.53,
                           0.68, 0.87, 1.12, 1.44, 1.85, 2.38, 3.06, 3.93, 5.05, 6.49,
                           8.34, 10.72, 13.78, 17.71, 22.76])

def read_last_flux(file_path):
    """Reads the last time-step row of flux values from file."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    data_lines = [line.strip() for line in lines if line.strip() and not line.startswith('#')]
    last_values = data_lines[-1].split()
    return np.array(list(map(float, last_values[1:])))  # skip time column

# Storage for plotting
eta_results = {}

# Compute deposition efficiencies
for label, folder in folders.items():
    inlet_flux = read_last_flux(f"{folder}/patch.inlet.dat")
    wall_flux = read_last_flux(f"{folder}/patch.walls.dat")
    eta = -wall_flux / inlet_flux
    eta_results[label] = eta

# Extract eta arrays
eta_shape_1 = eta_results['Shape factor = 1']
eta_shape_1_5 = eta_results['Shape factor = 1.5']

# Compute ratio of efficiencies
eta_ratio = np.divide(eta_shape_1,eta_shape_1_5, out=np.full_like(eta_shape_1, np.nan), where=eta_shape_1 != 0)

# Save η and ratio data
with open("Eta_ratio_vs_size.txt", "w") as f:
    f.write("Particle Size (µm)   η_Shape_1      η_Shape_1.5    η_ratio_1.5_to_1\n")
    f.write("-------------------   ------------    ------------    -----------------\n")
    for d, eta1, eta15, ratio in zip(particle_sizes, eta_shape_1, eta_shape_1_5, eta_ratio):
        f.write(f"{d:>17.5f}   {eta1:>12.6f}    {eta15:>12.6f}    {ratio:>17.6f}\n")

# Plotting
plt.figure(figsize=(8, 5), dpi=300)
for label in folders:
    plt.loglog(particle_sizes, eta_results[label], marker='o', label=label)

# Font settings
font_size = 10
tick_size = 10

plt.xlabel('Particle Size (μm)', fontsize=font_size, fontweight='bold')
plt.ylabel('Deposition Efficiency (η)', fontsize=font_size, fontweight='bold')

plt.xticks(fontsize=tick_size, fontweight='bold')
plt.yticks(fontsize=tick_size, fontweight='bold')

plt.legend(fontsize=10)

# Add simulation parameters as text on plot
plt.text(0.07, 0.8, 'D = 10 mm\nR$^*$ = 5.7\nRe = 1000\nDe = 419',
         fontsize=9, fontweight='bold',
         transform=plt.gca().transAxes, verticalalignment='top')

plt.tight_layout()

# Save the figure
plt.savefig('deposition_efficiency_plot.png', dpi=600)
plt.show()

