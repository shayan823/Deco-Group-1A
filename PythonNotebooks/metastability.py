from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from network import solve_ode_network, signal_processing_on_hopf, wc_ode_network
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
num_steps = 300000
# G = 0.5

mat = loadmat('../references/AAL78/C78.mat')
C = mat['C']
alpha = 0.2
C = alpha * C / np.max(C)

def metastability_vs_G():
    G_values = np.linspace(0, 1, 10)
    betas = [0.5, 1.5, 3]
    
    plt.figure(figsize=(10, 5))

    # Loop through the provided beta values
    for beta in betas:
        metastabilities = []

        # Calculate metastability for varying G
        for G in G_values:
            x_values, _ = solve_ode_network(num_steps, dt, a, omega, beta, C, G)
            _, _, analytical_x, _, _, _ = signal_processing_on_hopf(x_values, C, dt, f, G)
            _, metastability = calculate_kuramoto_and_metastability(analytical_x)
            metastabilities.append(metastability)

        # Plot metastability over G values for the current beta
        plt.plot(G_values, metastabilities, marker='o', label=f'beta = {beta}')

    # Finalize the plot with title, labels, and legend
    plt.title('Metastability for Stuart-Landau with varying G')
    plt.xlabel('G (Coupling Strength)')
    plt.ylabel('Metastability (Standard Deviation of R(t))')
    plt.legend()
    plt.tight_layout()

    # Save the combined plot
    plt.savefig(f"../Plots/python_plots/metstability/ms_vs_G_stuart_landau_tmax={dt*num_steps}_f={f}.png")

# def metastability_vs_G(beta):
#     G_values = np.linspace(0, 1, 10)
#     metastabilities = []

#     # Calculate metastability for varying G and plot the metastability
#     for G in G_values:
#         x_values, _ = solve_ode_network(num_steps, dt, a, omega, beta, C, G)
#         _, _, analytical_x, _, _, _ = signal_processing_on_hopf(x_values, C, dt, f, G)
#         _, metastability = calculate_kuramoto_and_metastability(analytical_x)
#         metastabilities.append(metastability)

#     # Plot metastability over G values
#     plt.figure(figsize=(10, 5))
#     plt.plot(G_values, metastabilities, marker='o')
#     plt.title('')
#     plt.xlabel('G (Coupling Strength)')
#     plt.ylabel('Metastability (Standard Deviation of R(t))')
#     plt.tight_layout()
#     plt.savefig(f"../Plots/python_plots/metstability/ms_vs_G_beta={beta}_tmax={dt*num_steps}_f={f}.png")

def metastability_vs_G_for_wc():
    p_i = -1
    G_values = np.linspace(0, 1, 10)
    n=3.5
    betas = [0.05,0.1,0.5]
    
    # Plot metastability over G values
    plt.figure(figsize=(10, 5))
    
    for beta in betas:
        metastabilities = []

        # Calculate metastability for varying G and plot the metastability
        for G in G_values:
            e_values, _ = wc_ode_network(p_i, num_steps, dt, beta, C, n, G)
            _, _, analytical_x, _, _, _ = signal_processing_on_hopf(e_values, C, dt, f, G)
            _, metastability = calculate_kuramoto_and_metastability(analytical_x)
            metastabilities.append(metastability)

        plt.plot(G_values, metastabilities, marker='o', label=f'beta = {beta}')

    plt.title('Metastability for Varying Coupling Strength G')
    plt.xlabel('G (Coupling Strength)')
    plt.ylabel('Metastability (Standard Deviation of R(t))')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../Plots/python_plots/wc_plots/metstability/ms_vs_G_wilson_cowan_tmax={dt*num_steps}_f={f}.png")


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

def metastability_vs_f_for_wc(beta):
    # Vary the frequencies from 4 to 28 with 12 points
    G=0.5
    f_values = [8,10,12,14,16,18,20]
    metastabilities_e = []
    metastabilities_i = []
    p_i = -1
    n = 3.5

    e_values, i_values = wc_ode_network(p_i, num_steps, dt, beta, C, n, G)
    # Calculate metastability for varying G and plot the metastability
    for f in f_values:
        _, _, analytical_e, _, _, _ = signal_processing_on_hopf(e_values, C, dt, f, G)
        _, metastability_e = calculate_kuramoto_and_metastability(analytical_e)
        metastabilities_e.append(metastability_e)

    # Plot metastability over G values
    plt.figure(figsize=(10, 5))
    plt.plot(f_values, metastabilities_e, marker='o', color='red', label='Excitatory signal')
    plt.title('Metastability for Varying Carrier frequency f')
    plt.xlabel('f (Carrier Frequency)')
    plt.ylabel('Metastability (Standard Deviation of R(t))')
    plt.legend()
    plt.tight_layout()
    # plt.savefig(f"../Plots/python_plots/metstability/ms_vs_f_beta={beta}_tmax={dt*num_steps}.png")
    plt.savefig(f"../Plots/python_plots/wc_plots/metstability/ms_vs_f_beta={beta}_tmax={dt*num_steps}.png")

# metastability_vs_G()
# metastability_vs_G(3)
# metastability_vs_f(0.5)
# metastability_vs_f(3)

metastability_vs_G_for_wc()
# metastability_vs_G_for_wc(0.5)
# metastability_vs_f_for_wc(0.1)
# metastability_vs_f_for_wc(0.5)

