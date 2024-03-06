from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import hilbert, butter, filtfilt
from scipy.signal import hilbert, hilbert2, savgol_filter, find_peaks
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
from scipy.signal import iirdesign, sosfilt
from tqdm import tqdm
from scipy.signal import butter, lfilter

mat = loadmat('../references/AAL78/C78.mat')
C = mat['C']
alpha = 0.2
C = alpha * C / np.max(C)


def solve_ode_network(num_steps, dt, a, omega, beta, C, G=0.5, seed=42):
    print("Solving ODE Network...")
    # Initialize the arrays to store the x, y values for all neurons
    np.random.seed(seed)

    total_nodes = C.shape[0]

    x_values = np.zeros((num_steps, total_nodes))
    y_values = np.zeros((num_steps, total_nodes))
    
    # Initialize x and y vectors with zeros
    x = np.ones(total_nodes)*0.5
    y = np.ones(total_nodes)*0.5

    # Pre-generate noise for all neurons and all time steps
    noise = np.random.randn(num_steps, 78) * np.sqrt(dt)

    for step in range(num_steps):
        # Generate a single random noise term for each neuron
        # Extract the noise for the current step
        current_noise = noise[step, :]


        # Calculate the dxdt and dydt using the equations provided
        dxdt = (a - x**2 - y**2)*x - omega*y + G * np.dot(C, (x - x[:, None]).T).diagonal() 
        dydt = (a - x**2 - y**2)*y + omega*x + G * np.dot(C, (y - y[:, None]).T).diagonal() 
        
        # Update the x and y values
        x += dxdt * dt + beta * current_noise
        y += dydt * dt + beta * current_noise
        
        # Store the values
        x_values[step, :] = x
        y_values[step, :] = y

    
    return x_values, y_values

def wc_ode_network(p_i, num_steps, dt, beta, C, n, G=0.5, seed=42):

    tau_e = 0.06/n
    tau_i = 0.16/n
    w_ee = 10
    w_ei = 12
    w_ii = 10
    w_ie = 10
    p_e = 0.5122

    total_nodes = C.shape[0]

    np.random.seed(seed)
    e_values = np.zeros((num_steps, total_nodes))
    i_values = np.zeros((num_steps, total_nodes))

    e = np.ones(total_nodes)*1
    i = np.ones(total_nodes)*1

    f_a = lambda x: 1 / (1 + np.exp(-x))

    # Pre-generate noise for all neurons and all time steps
    noise = np.random.randn(num_steps, total_nodes) * np.sqrt(dt)
    
    for step in range(num_steps):

        current_noise = noise[step, :]

        coupling = G*(np.sum(C, axis=1))

        dedt = (-e + f_a(w_ee*e - w_ei*i + p_e + coupling))/tau_e
        didt = (-i + f_a(w_ie*e - w_ii*i + p_i))/tau_i

        e += dedt * dt + beta * current_noise
        i += didt * dt + beta * current_noise

        e_values[step, :] = e
        i_values[step, :] = i
    return e_values, i_values

# Function to apply band-pass filter
def bandpass_filter(data, lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    y = lfilter(b, a, data)
    return y

# Function to apply low-pass filter
def lowpass_filter(data, cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = lfilter(b, a, data)
    return y

def bandpass_filter_auto_order(data, lowcut, highcut, fs, gpass=3, gstop=40):
    # Normalized passband and stopband edges
    wp = [lowcut / (0.5 * fs), highcut / (0.5 * fs)]
    ws = [(lowcut - lowcut*0.1) / (0.5 * fs), (highcut + highcut*0.1) / (0.5 * fs)]
    
    # Design the filter using iirdesign
    sos = iirdesign(wp, ws, gpass, gstop, analog=False, ftype='butter', output='sos')
    
    # Apply the filter
    y = sosfilt(sos, data)
    return y

def butter_bandpass(lowcut, highcut, fs, order=5):
    return butter(order, [lowcut, highcut], fs=fs, btype='band')

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def smoothening_envelope(amplitude_env):
    # Initialize the smooth outline matrix with the same shape as amplitude_env
    smooth_outline = np.zeros_like(amplitude_env)
    
    # Iterate over each node and compute the smooth outline
    for i in range(amplitude_env.shape[1]):
        # Find peaks for the current node's amplitude envelope
        peaks, _ = find_peaks(amplitude_env[:, i])
        
        # Handle the case when no peaks are found
        if len(peaks) == 0:
            # If no peaks are found, we could choose to return the original signal,
            # or handle this case according to the specific requirements of the analysis
            smooth_outline[:, i] = amplitude_env[:, i]
            continue
        x_values = np.arange(amplitude_env.shape[0])
        # Interpolate between peaks to create a smoother outline
        peak_interpolation = interp1d(peaks, amplitude_env[peaks, i], kind='cubic', bounds_error=False, fill_value='extrapolate')
        
        # Use the interpolation function to generate the smooth outline for the current node
        smooth_outline[:, i] = peak_interpolation(x_values)
    return smooth_outline

def signal_processing_on_hopf(signal, C, dt, f=12, G=0.5):   
    print(f"Processing signal for frequency {f}Hz...")
    # Filter parameters
    fs = 1/dt
    band_lowcut = f - 2
    band_highcut = f + 2
    lowpass_cutoff = 0.2

    # Preallocate arrays for filtered signals
    filtered_signal = np.zeros_like(signal)

    # Band-pass filter for each region
    for i in range(78):
        # filtered_signal[:, i] = bandpass_filter(signal[:, i], band_lowcut, band_highcut, fs)
        # filtered_signal[:, i] = bandpass_filter_auto_order(signal[:, i], band_lowcut, band_highcut, fs)
        filtered_signal[:, i] = butter_bandpass_filter(signal[:, i], band_lowcut, band_highcut, fs)
        
    # Calculate the Hilbert transform to get the amplitude envelope
    amplitude_env = np.abs(hilbert(filtered_signal))

    # Smoothening the envelope
    smoothened_env = smoothening_envelope(filtered_signal)

    # Low-pass filter the amplitude envelope
    analytical_signal = hilbert(smoothened_env)
    ultra_slow = np.zeros_like(amplitude_env) #smoothened_env

    for i in range(78):
        # ultra_slow[:, i] = lowpass_filter(smoothened_env[:, i], lowpass_cutoff, fs) # smoothened_env
        ultra_slow[:, i] = butter_lowpass_filter(amplitude_env[:, i], lowpass_cutoff, fs) # smoothened_env

    return signal, filtered_signal, analytical_signal, amplitude_env, smoothened_env, ultra_slow


def compute_envelope_fc(ultra_slow_x):
    print("Computing Envelope Functional Connectivity...")
    # Number of regions
    n_regions = ultra_slow_x.shape[1]
    
    # Initialize the FC matrix
    envelope_fc = np.zeros((n_regions, n_regions))
    
    # Compute Pearson's correlation for each pair of brain regions
    for i in range(n_regions):
        for j in range(n_regions):
            if i != j:
                corr, _ = pearsonr(ultra_slow_x[:, i], ultra_slow_x[:, j])
                envelope_fc[i, j] = corr
            else:
                envelope_fc[i, j] = 1  # Self-correlation is always 1
                
    return envelope_fc