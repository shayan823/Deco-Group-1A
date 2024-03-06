# Import necessary libraries
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from network import solve_ode_network, signal_processing_on_hopf, smoothening_envelope, compute_envelope_fc
from scipy.stats import pearsonr
from scipy.signal import hilbert, butter, filtfilt, hilbert2, savgol_filter, find_peaks
from scipy.interpolate import interp1d
from matplotlib.ticker import FixedLocator, FixedFormatter

# Load connectivity matrix
mat = loadmat('../references/AAL78/C78.mat')
C = mat['C']

# Scale the connectivity matrix
alpha = 0.2
C = alpha * C / np.max(C)

num_steps = 25000  # total steps
dt = 0.001  # time step
a = 0.0  # intrinsic dynamics parameters for each neuron

def get_v_t_data(beta, f):
    G=0.5
    omega = 2 * np.pi * 12  # angular frequencies for each neuron

    x_values, _ = solve_ode_network(num_steps, dt, a, omega, beta, C, G)

    signal = x_values # r_values
    _, _, analytical_x, _, _, _ =  signal_processing_on_hopf(signal, C, dt, f=12, G=0.5)

    instantaneous_phase = np.angle(analytical_x)
    print(instantaneous_phase.shape)
    # Number of time points and brain regions
    num_time_points, num_brain_regions = instantaneous_phase.shape

    # Step 1: Calculate V(t) for each time step
    V_t = []
    for t in range(num_time_points):
        phase_diff = np.abs(instantaneous_phase[t, :, None] - instantaneous_phase[t, None, :])
        cos_phase_diff = np.cos(phase_diff)
        # We exclude the diagonal to not include the coherence of a node with itself
        np.fill_diagonal(cos_phase_diff, 0)
        V_t.append(cos_phase_diff[np.triu_indices(num_brain_regions, k=1)])
    V_t = np.array(V_t)
    print(V_t.shape)
    np.save(f'./data/v_t_beta={beta}_tmax={dt*num_steps}_f={f}.npy', V_t)

def calculate_ccd(beta, f):
    V_t = np.load(f'./data/v_t_beta={beta}_tmax={dt*num_steps}_f={f}.npy')
    # Step 2: Calculate the CCD matrix
    CCD = np.corrcoef(V_t)
    print(CCD.shape)
    np.save(f'./data/ccd_beta={beta}_tmax={dt*num_steps}_f={f}.npy', CCD)

def plot_ccd_matrix(beta, f, dt, num_steps):
    # Load the CCD data
    CCD = np.load(f'./data/ccd_beta={beta}_tmax={dt*num_steps}_f={f}.npy')
    print(f"CCD Matrix shape: {CCD.shape}")

    # Plotting the CCD matrix
    plt.figure(figsize=(6, 6))  # Adjust the figure size as needed
    im = plt.imshow(CCD, cmap='turbo', origin='lower', aspect='auto')
    plt.colorbar(im, label='Cosine Similarity')
    plt.title('CCD Matrix')
    plt.xlabel('Time (t2)')
    plt.ylabel('Time (t1)')
    plt.gca().invert_yaxis()  # Invert the y-axis

    # Customize tick labels to be divided by 1000
    xticks = plt.gca().get_xticks()
    yticks = plt.gca().get_yticks()
    plt.gca().xaxis.set_major_locator(FixedLocator(xticks))
    plt.gca().yaxis.set_major_locator(FixedLocator(yticks))
    plt.gca().set_xticklabels([f"{x/1000:.0f}" for x in xticks])
    plt.gca().set_yticklabels([f"{y/1000:.0f}" for y in yticks])

    plt.tight_layout()
    plt.savefig(f"../Plots/python_plots/ccd/beta={beta}_tmax={dt*num_steps}_f={f}.png")

def plot_ccd_distribution(beta, f, dt, num_steps):
    # Load the CCD data
    CCD = np.load(f'./data/ccd_beta={beta}_tmax={dt*num_steps}_f={f}.npy')
    CCD_values = CCD[np.triu_indices_from(CCD, k=1)]

    # Plotting the distribution of CCD values
    plt.figure(figsize=(6, 4))  # Adjust the figure size as needed
    plt.hist(CCD_values, bins=100, color='blue', alpha=0.7)
    plt.title(f'Distribution of CCD Values for beta={beta}')
    plt.xlabel('Cosine Similarity')
    plt.ylabel('Count')

    plt.tight_layout()
    plt.savefig(f"../Plots/python_plots/ccd/ccd_distbeta={beta}_tmax={dt*num_steps}_f={f}.png")

# get_v_t_data(0.5, 12)
# calculate_ccd(0.5, 12)
# plot_ccd_matrix(0.5, 12, dt, num_steps)
# plot_ccd_distribution(0.5, 12, dt, num_steps)

get_v_t_data(3, 12)
calculate_ccd(3, 12) 
plot_ccd_matrix(3, 12, dt, num_steps)
plot_ccd_distribution(3, 12, dt, num_steps)