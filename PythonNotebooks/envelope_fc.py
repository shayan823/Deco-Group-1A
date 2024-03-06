"""
This script loads a connectivity matrix from a .mat file and processes it to simulate and analyze the
functional connectivity based on the Hopf model and signal processing techniques. It generates envelope
functional connectivity (FC) matrices for a range of frequencies and plots the results as heatmaps.
"""

# Import necessary libraries
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from network import solve_ode_network, signal_processing_on_hopf, compute_envelope_fc
from scipy.stats import pearsonr
from scipy.signal import hilbert, butter, filtfilt, hilbert2, savgol_filter, find_peaks
from scipy.interpolate import interp1d

# Load connectivity matrix
mat = loadmat('../references/AAL78/C78.mat')
C = mat['C']

# Scale the connectivity matrix
alpha = 0.2
C = alpha * C / np.max(C)

# Simulation parameters
G = 0.5
num_steps = 300000  # total steps
dt = 0.001  # time step
a = 0.0  # intrinsic dynamics parameters for each neuron
# beta = 0.5  # noise level

# frequencies = np.arange(12, 15, 4)
frequencies = [8, 12, 16]

def multiple_env_fc(beta):
    # Set up the plot grid
    num_frequencies = len(frequencies)
    cols = 3  # Define number of columns for subplots
    rows = num_frequencies // cols + (num_frequencies % cols > 0)  # Calculate required number of rows
    plt.figure(figsize=(cols * 5, rows * 4))  # Adjust the size as needed

    # Determine the global min and max correlation values for consistent colormap scaling
    # min_corr = np.inf
    # max_corr = -np.inf
    # for f in frequencies:
    #     omega = 2 * np.pi * f  # angular frequencies for each neuron
    #     x_values, y_values = solve_ode_network(num_steps, dt, a, omega, beta, C, G)
    #     _, _, _, _, _, ultra_slow_x = signal_processing_on_hopf(x_values, C, dt, f, G)
    #     envelope_fc = compute_envelope_fc(ultra_slow_x)
    #     min_corr = min(min_corr, envelope_fc.min())
    #     max_corr = max(max_corr, envelope_fc.max())

    # Generate and plot Envelope FC matrices for each frequency
    for i, f in enumerate(frequencies):
        omega = 2 * np.pi * f  # angular frequencies for each neuron
        x_values, _ = solve_ode_network(num_steps, dt, a, omega, beta, C, G)
        
        _, _, _, _, _, ultra_slow_x = signal_processing_on_hopf(x_values, C, dt, f, G)
        envelope_fc = compute_envelope_fc(ultra_slow_x)

        # Plotting
        plt.subplot(rows, cols, i + 1)
        sns.heatmap(envelope_fc, cmap='turbo', square=True, cbar_kws={'label': 'Pearson Correlation'})
        plt.title(f'{f-2}-{f+2} Hz')
        plt.xlabel('Brain Regions')
        plt.ylabel('Brain Regions')

    # Save the figure
    plt.tight_layout()
    plt.savefig(f"../Plots/python_plots/envelope_fc/multi_fc/multiple_fc_beta={beta}_tmax={dt*num_steps}_v1.png")


def single_freq_envelope_fc(beta):
    print("Starting single_freq_envelope_fc function...")
    omega = 2 * np.pi * 12  # angular frequencies for each neuron
    print("Calling solve_ode_network...")
    x_values, y_values = solve_ode_network(num_steps, dt, a, omega, beta, C, G)

    print("Setting up plot grid...")
    num_frequencies = len(frequencies)
    cols = 3  # Define number of columns for subplots
    rows = num_frequencies // cols + (num_frequencies % cols > 0)  # Calculate required number of rows

    plt.figure(figsize=(cols * 5, rows * 4))  # Adjust the size as needed

    min_corr = np.inf
    max_corr = -np.inf

    for f in frequencies:
        print(f"Processing frequency: {f}Hz")
        _, _, _, _, _, ultra_slow_x = signal_processing_on_hopf(x_values, C, dt, f, G)
        envelope_fc = compute_envelope_fc(ultra_slow_x)
        min_corr = min(min_corr, envelope_fc.min())
        max_corr = max(max_corr, envelope_fc.max())

    for i, f in enumerate(frequencies):
        print(f"Plotting Envelope FC Matrix for {f}Hz...")
        _, _, _, _, _, ultra_slow_x = signal_processing_on_hopf(x_values, C, dt, f, G)
        envelope_fc = compute_envelope_fc(ultra_slow_x)

        plt.subplot(rows, cols, i + 1)
        sns.heatmap(envelope_fc, cmap='turbo', square=True, cbar_kws={'label': 'Pearson Correlation'}, vmin=min_corr, vmax=max_corr)
        plt.title(f'{f-2}-{f+2} Hz')
        plt.xlabel('Brain Regions')
        plt.ylabel('Brain Regions')

    plt.tight_layout()
    plt.savefig(f"../Plots/python_plots/envelope_fc/beta={beta}_tmax={dt*num_steps}_tmp.png")
    print("Plot saved.")


def plot_single_for_12(beta):
    f=12
    omega = 2 * np.pi * 12
    x_values, _ = solve_ode_network(num_steps, dt, a, omega, beta, C, G)
    # Save the simulation results
    np.save(f'./data/x_values_{f}Hz_{beta}.npy', x_values)

    print("Signal Processing")

    cols = 1  # Define number of columns for subplots
    rows = 1 // cols + (1 % cols > 0)
    plt.figure(figsize=(cols * 6, rows * 6))


    # Load the simulation results
    x_values = np.load(f'./data/x_values_{f}Hz_{beta}.npy')
    
    _, _, _, _, _, ultra_slow_x = signal_processing_on_hopf(x_values, C, dt, f, G)
    envelope_fc = compute_envelope_fc(ultra_slow_x)

    # Plotting
    sns.heatmap(envelope_fc, cmap='turbo', square=True, cbar_kws={'label': 'Pearson Correlation'})
    plt.title(f"Envelope FC Matrix - {f} Hz")
    plt.xlabel('Brain Regions')
    plt.ylabel('Brain Regions')

    plt.tight_layout()
    plt.savefig(f"../Plots/python_plots/envelope_fc/beta={beta}_tmax={dt*num_steps}_f={f}.png")



# single_freq_envelope_fc(0.5)
# single_freq_envelope_fc(3)
multiple_env_fc(0.5)
multiple_env_fc(3)
# plot_single_for_12(0.5)
# plot_single_for_12(3)