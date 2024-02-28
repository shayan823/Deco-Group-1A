from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import hilbert, butter, filtfilt
from scipy.signal import hilbert, hilbert2, savgol_filter, find_peaks
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
from scipy.signal import iirdesign, sosfilt

mat = loadmat('../references/AAL78/C78.mat')
C = mat['C']

def solve_ode_network(num_steps, dt, a, omega, beta, C, G=0.5, seed=42):
    # Initialize the arrays to store the x, y values for all neurons
    np.random.seed(seed)
    x_values = np.zeros((num_steps, 78))
    y_values = np.zeros((num_steps, 78))
    
    # Initialize x and y vectors with zeros
    x = np.ones(78)*0.5
    y = np.ones(78)*0.5

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

# Function to apply band-pass filter
def bandpass_filter(data, lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y

# Function to apply low-pass filter
def lowpass_filter(data, cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
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

# def smoothening_envelope(amplitude_env, smoothing_window_length=51, polyorder=3, prominence=None, width=None, distance=None):
#     smooth_outline = np.zeros_like(amplitude_env)
    
#     # Apply a Savitzky-Golay filter to the amplitude envelope for pre-smoothing
#     for i in range(amplitude_env.shape[1]):
#         smoothed_env = savgol_filter(amplitude_env[:, i], smoothing_window_length, polyorder)

#         # Find peaks in the smoothed amplitude envelope
#         peaks, _ = find_peaks(smoothed_env, prominence=prominence, width=width, distance=distance)

#         # If no peaks are found, return the original signal or handle as required
#         if len(peaks) == 0:
#             smooth_outline[:, i] = amplitude_env[:, i]
#             continue

#         # Create the x-values for interpolation
#         x_values = np.arange(amplitude_env.shape[0])

#         # Interpolate between peaks to create a smoother outline
#         # If cubic interpolation fails due to insufficient points, fall back to linear
#         try:
#             peak_interpolation = interp1d(peaks, smoothed_env[peaks], kind='cubic', bounds_error=False, fill_value='extrapolate')
#         except ValueError:
#             peak_interpolation = interp1d(peaks, smoothed_env[peaks], kind='linear', bounds_error=False, fill_value='extrapolate')
        
#         # Generate the smooth outline for the current node
#         smooth_outline[:, i] = peak_interpolation(x_values)

    return smooth_outline

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
    # Filter parameters
    fs = 1/dt
    band_lowcut = f - 2
    band_highcut = f + 2
    lowpass_cutoff = 0.2

    # Preallocate arrays for filtered signals
    filtered_signal = np.zeros_like(signal)

    # Band-pass filter for each region
    for i in range(78):
        filtered_signal[:, i] = bandpass_filter(signal[:, i], band_lowcut, band_highcut, fs)
        # filtered_signal[:, i] = bandpass_filter_auto_order(signal[:, i], band_lowcut, band_highcut, fs)
        
    # Calculate the Hilbert transform to get the amplitude envelope
    amplitude_env = np.abs(hilbert(filtered_signal))

    # Smoothening the envelope
    smoothened_env = smoothening_envelope(amplitude_env)

    # Low-pass filter the amplitude envelope
    analytical_signal = hilbert(smoothened_env)
    ultra_slow = np.zeros_like(amplitude_env)

    for i in range(78):
        # ultra_slow[:, i] = lowpass_filter(amplitude_env[:, i], lowpass_cutoff, fs)
        ultra_slow[:, i] = lowpass_filter(smoothened_env[:, i], lowpass_cutoff, fs)

    return signal, filtered_signal, analytical_signal, amplitude_env, smoothened_env, ultra_slow

def compute_envelope_fc(ultra_slow_x):
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