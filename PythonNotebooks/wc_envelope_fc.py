from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from network import wc_ode_network, signal_processing_on_hopf, compute_envelope_fc
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
p_i = -1

frequencies = [12, 16, 20]

def single_freq_envelope_fc(beta):
    print("Starting single_freq_envelope_fc function...")
    omega = 2 * np.pi * 12  # angular frequencies for each neuron
    print("Calling solve_ode_network...")
    # x_values, y_values = solve_ode_network(num_steps, dt, a, omega, beta, C, G)
    e_values, _ = wc_ode_network(p_i, num_steps, dt, beta, C, n=3.5, G=0.5)
    # _, i_values = wc_ode_network(p_i, num_steps, dt, beta, C, n=3.5, G=0.5)
    signal = e_values

    print("Setting up plot grid...")
    num_frequencies = len(frequencies)
    cols = 3  # Define number of columns for subplots
    rows = num_frequencies // cols + (num_frequencies % cols > 0)  # Calculate required number of rows

    plt.figure(figsize=(cols * 5, rows * 4))  # Adjust the size as needed

    # min_corr = np.inf
    # max_corr = -np.inf

    # for f in frequencies:
    #     print(f"Processing frequency: {f}Hz")
    #     _, _, _, _, _, ultra_slow_x = signal_processing_on_hopf(e_values, C, dt, f, G)
    #     envelope_fc = compute_envelope_fc(ultra_slow_x)
    #     min_corr = min(min_corr, envelope_fc.min())
    #     max_corr = max(max_corr, envelope_fc.max())

    for i, f in enumerate(frequencies):
        print(f"Plotting Envelope FC Matrix for {f}Hz...")
        _, _, _, _, _, ultra_slow_x = signal_processing_on_hopf(signal, C, dt, f, G)
        envelope_fc = compute_envelope_fc(ultra_slow_x)

        plt.subplot(rows, cols, i + 1)
        # sns.heatmap(envelope_fc, cmap='turbo', square=True, cbar_kws={'label': 'Pearson Correlation'}, vmin=min_corr, vmax=max_corr)
        sns.heatmap(envelope_fc, cmap='turbo', square=True, cbar_kws={'label': 'Pearson Correlation'})
        plt.title(f'{f-2}-{f+2} Hz')
        plt.xlabel('Brain Regions')
        plt.ylabel('Brain Regions')

    plt.tight_layout()
    plt.savefig(f"../Plots/wc_plots/envelope_fc/e_signal_beta={beta}_tmax={dt*num_steps}.png")
    # plt.savefig(f"../Plots/wc_plots/envelope_fc/i_signal_beta={beta}_tmax={dt*num_steps}.png")
    print("Plot saved.")

single_freq_envelope_fc(0.05)
# single_freq_envelope_fc(3)