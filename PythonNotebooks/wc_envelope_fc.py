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
num_steps = 400000  # total steps
dt = 0.001  # time step
p_i = -1

frequencies = [12, 16, 20]

# def single_freq_envelope_fc(beta):
#     print("Starting single_freq_envelope_fc function...")
#     omega = 2 * np.pi * 12  # angular frequencies for each neuron
#     print("Calling solve_ode_network...")
#     # x_values, y_values = solve_ode_network(num_steps, dt, a, omega, beta, C, G)
#     e_values, _ = wc_ode_network(p_i, num_steps, dt, beta, C, n=3.5, G=0.5)
#     # _, i_values = wc_ode_network(p_i, num_steps, dt, beta, C, n=3.5, G=0.5)
#     signal = e_values

#     print("Setting up plot grid...")
#     num_frequencies = len(frequencies)
#     cols = 3  # Define number of columns for subplots
#     rows = num_frequencies // cols + (num_frequencies % cols > 0)  # Calculate required number of rows

#     plt.figure(figsize=(cols * 5, rows * 4))  # Adjust the size as needed

#     # min_corr = np.inf
#     # max_corr = -np.inf

#     # for f in frequencies:
#     #     print(f"Processing frequency: {f}Hz")
#     #     _, _, _, _, _, ultra_slow_x = signal_processing_on_hopf(e_values, C, dt, f, G)
#     #     envelope_fc = compute_envelope_fc(ultra_slow_x)
#     #     min_corr = min(min_corr, envelope_fc.min())
#     #     max_corr = max(max_corr, envelope_fc.max())

#     for i, f in enumerate(frequencies):
#         print(f"Plotting Envelope FC Matrix for {f}Hz...")
#         _, _, _, _, _, ultra_slow_x = signal_processing_on_hopf(signal, C, dt, f, G)
#         envelope_fc = compute_envelope_fc(ultra_slow_x)

#         plt.subplot(rows, cols, i + 1)
#         # sns.heatmap(envelope_fc, cmap='turbo', square=True, cbar_kws={'label': 'Pearson Correlation'}, vmin=min_corr, vmax=max_corr)
#         sns.heatmap(envelope_fc, cmap='turbo', square=True, cbar_kws={'label': 'Pearson Correlation'})
#         plt.title(f'{f-2}-{f+2} Hz')
#         plt.xlabel('Brain Regions')
#         plt.ylabel('Brain Regions')

#     plt.tight_layout()
#     plt.savefig(f"../Plots/wc_plots/envelope_fc/e_signal_beta={beta}_tmax={dt*num_steps}.png")
#     # plt.savefig(f"../Plots/wc_plots/envelope_fc/i_signal_beta={beta}_tmax={dt*num_steps}.png")
#     print("Plot saved.")

# single_freq_envelope_fc(0.05)
# single_freq_envelope_fc(3)


# def calculate_mean_correlation(A, B):
#     """
#     Calculate mean Pearson correlation between two matrices.
#     """
#     # Flatten the matrices
#     A_flat = A.flatten()
#     B_flat = B.flatten()
#     # Calculate Pearson correlation coefficient
#     correlation, _ = pearsonr(A_flat, B_flat)
#     return correlation

def plot_envelope_fc(frequencies, betas, C, dt, num_steps, p_i):
    cols = len(frequencies)  # Columns for different frequencies
    rows = len(betas)  # Rows for different beta values

    # Creating figures for excitatory and inhibitory signals
    fig_e, axs_e = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4), sharex=True, sharey=True)
    fig_i, axs_i = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4), sharex=True, sharey=True)

    # Loop through betas and frequencies to fill the subplots
    for row, beta in enumerate(betas):
        for col, f in enumerate(frequencies):
            print(f"Processing beta={beta}, frequency={f}Hz...")

            e_values, i_values = wc_ode_network(p_i, num_steps, dt, beta, C, n=3.5, G=0.5)
            ultra_slow_e, envelope_fc_e = signal_processing_and_fc(e_values, C, dt, f)
            ultra_slow_i, envelope_fc_i = signal_processing_and_fc(i_values, C, dt, f)

            # Plot for excitatory values
            ax = axs_e[row, col] if rows > 1 else axs_e[col]
            sns.heatmap(envelope_fc_e, ax=ax, cmap='turbo', square=True, cbar_kws={'label': 'Pearson Correlation'})
            if col == 0:
                ax.set_ylabel(f'Beta={beta}')
            if row == 0:
                ax.set_title(f'{f-2}-{f+2} Hz')

            # Plot for inhibitory values
            ax = axs_i[row, col] if rows > 1 else axs_i[col]
            sns.heatmap(envelope_fc_i, ax=ax, cmap='turbo', square=True, cbar_kws={'label': 'Pearson Correlation'})
            if col == 0:
                ax.set_ylabel(f'Beta={beta}')
            if row == 0:
                ax.set_title(f'{f-2}-{f+2} Hz')

    fig_e.tight_layout()
    fig_i.tight_layout()
    fig_e.suptitle('Excitatory Envelope FC', y=1.02)
    fig_i.suptitle('Inhibitory Envelope FC', y=1.02)
    fig_e.savefig(f"../Plots/python_plots/wc_plots/envelope_fc/envelope_fc_combined_excitatory_{dt*num_steps}.png")
    fig_i.savefig(f"../Plots/python_plots/wc_plots/envelope_fc/envelope_fc_combined_inhibitory{dt*num_steps}.png")

    # plot_mean_correlation_scores(e_corrs, 'e_values', C)
    # plot_mean_correlation_scores(i_corrs, 'i_values', C)

# def plot_mean_correlation_scores(corrs, label, C):
#     fig, ax = plt.subplots(figsize=(10, 5))
#     for beta, f, corr in corrs:
#         ax.plot(f, corr, 'o-', label=f'beta={beta}')
#     ax.set_xlabel('Frequency (Hz)')
#     ax.set_ylabel('Mean Correlation Score')
#     ax.set_title(f'Mean Correlation Scores with C for {label}')
#     ax.legend()
#     plt.tight_layout()
#     plt.savefig(f"../Plots/python_plots/wc_plots/envelope_fc/mean_correlation_scores_{label}_with_C.png")


def signal_processing_and_fc(values, C, dt, f):
    # Placeholder for actual signal processing and FC calculation function
    _, _, _, _, _, ultra_slow_x = signal_processing_on_hopf(values, C, dt, f, G)
    envelope_fc = compute_envelope_fc(ultra_slow_x)
    return ultra_slow_x, envelope_fc

# Example usage
betas = [0.05, 0.1, 0.5]
p_i = 0.1  # Placeholder for initial condition parameter
plot_envelope_fc(frequencies, betas, C, dt, num_steps, p_i)