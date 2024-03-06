from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from network import solve_ode_network, signal_processing_on_hopf, smoothening_envelope, compute_envelope_fc
from utility import plot_x_values
from scipy.stats import pearsonr
from scipy.signal import hilbert, butter, filtfilt
from scipy.signal import hilbert, hilbert2, savgol_filter, find_peaks
from scipy.interpolate import interp1d

# Function to calculate Kuramoto order parameter R(t) and metastability
def calculate_kuramoto_and_metastability(analytical_x):
    # Calculate the instantaneous phase φk(t) for each brain region
    instantaneous_phase = np.angle(analytical_x)

    # Calculate the Kuramoto order parameter R(t) over time
    R_t = np.abs(np.sum(np.exp(1j * instantaneous_phase), axis=1))/78

    # Calculate metastability as the standard deviation of the Kuramoto order parameter
    metastability = np.std(R_t)
    
    return R_t, metastability

# beta = 3 # noise level
dt = 0.001
a = 0
f=12
omega = 2*np.pi*f
num_steps = 100000
# G = 0.5

mat = loadmat('../references/AAL78/C78.mat')
C = mat['C']
alpha = 0.2
C = alpha * C / np.max(C)


def metastability_vs_G(beta):
    G_values = np.linspace(0, 1, 10)
    metastabilities = []

    # Calculate metastability for varying G and plot the metastability
    for G in G_values:
        x_values, _ = solve_ode_network(num_steps, dt, a, omega, beta, C, G)
        _, _, analytical_x, _, _, _ = signal_processing_on_hopf(x_values, C, dt, f, G)
        _, metastability = calculate_kuramoto_and_metastability(analytical_x)
        metastabilities.append(metastability)

    # Plot metastability over G values
    plt.figure(figsize=(10, 5))
    plt.plot(G_values, metastabilities, marker='o')
    plt.title('Metastability for Varying Coupling Strength G')
    plt.xlabel('G (Coupling Strength)')
    plt.ylabel('Metastability (Standard Deviation of R(t))')
    plt.tight_layout()
    plt.savefig(f"../Plots/python_plots/metstability/ms_vs_G_beta={beta}_tmax={dt*num_steps}.png")

def metastability_vs_f(beta):
    # Vary the frequencies from 4 to 28 with 12 points
    G=0.5
    f_values = [8,10,12,14,16,18,20]
    metastabilities = []
    metastabilities_single = []

    # Calculate metastability for varying G and plot the metastability
    for f in f_values:
        omega = 2 * np.pi * 12
        x_values, _ = solve_ode_network(num_steps, dt, a, omega, beta, C, G)
        _, _, analytical_x, _, _, _ = signal_processing_on_hopf(x_values, C, dt, f, G)
        _, metastability = calculate_kuramoto_and_metastability(analytical_x)
        metastabilities_single.append(metastability)

        omega = 2 * np.pi * f
        x_values, _ = solve_ode_network(num_steps, dt, a, omega, beta, C, G)
        _, _, analytical_x, _, _, _ = signal_processing_on_hopf(x_values, C, dt, f, G)
        _, metastability = calculate_kuramoto_and_metastability(analytical_x)
        metastabilities.append(metastability)

    # Plot metastability over G values
    plt.figure(figsize=(10, 5))
    plt.plot(f_values, metastabilities_single, marker='o', color='red', label='Single frequency model')
    plt.plot(f_values, metastabilities, marker='o', color='blue', label='Multi frequency model')
    plt.title('Metastability for Varying Carrier frequency f')
    plt.xlabel('f (Carrier Frequency)')
    plt.ylabel('Metastability (Standard Deviation of R(t))')
    plt.legend()
    plt.tight_layout()
    # plt.savefig(f"../Plots/python_plots/metstability/ms_vs_f_beta={beta}_tmax={dt*num_steps}.png")
    plt.savefig(f"../Plots/python_plots/metstability/ms_vs_f_beta={beta}_tmax={dt*num_steps}.png")

# metastability_vs_G(0.5)
# metastability_vs_G(3)
metastability_vs_f(0.5)
metastability_vs_f(3)