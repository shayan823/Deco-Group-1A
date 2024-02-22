from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import hilbert, butter, filtfilt
from scipy.signal import hilbert, hilbert2, savgol_filter, find_peaks
from scipy.interpolate import interp1d
from scipy.stats import pearsonr

mat = loadmat('../references/AAL78/C78.mat')
C = mat['C']

def solve_ode_network(num_steps, dt, a, omega, beta, C, G=0.5):
    # Initialize the arrays to store the x, y values for all neurons
    x_values = np.zeros((num_steps, 78))
    y_values = np.zeros((num_steps, 78))
    
    # Initialize x and y vectors with zeros
    x = np.ones(78)*0.5
    y = np.ones(78)*0.5

    for step in range(num_steps):
        # Generate a single random noise term for each neuron
        noise = np.random.randn(78, 1) * np.sqrt(dt)
        noise = noise.reshape(78)

        # Calculate the dxdt and dydt using the equations provided
        dxdt = (a - x**2 - y**2)*x - omega*y + G * np.dot(C, (x - x[:, None]).T).diagonal() 
        dydt = (a - x**2 - y**2)*y + omega*x + G * np.dot(C, (y - y[:, None]).T).diagonal() 
        
        # Update the x and y values
        x += dxdt * dt + beta*noise
        y += dydt * dt + beta*noise
        
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

def smoothening_envelope(amplitude_env, band_pass):
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
        
        # Interpolate between peaks to create a smoother outline
        peak_interpolation = interp1d(peaks, amplitude_env[peaks, i], kind='linear', bounds_error=False, fill_value='extrapolate')
        
        # Use the interpolation function to generate the smooth outline for the current node
        smooth_outline[:, i] = peak_interpolation(np.arange(len(band_pass)))
    return smooth_outline

def signal_processing_on_hopf(C, f=12, G=0.5):   
    num_steps = 10000  # total steps
    dt = 0.001  # time step
    a = 0.0  # intrinsic dynamics parameters for each neuron
    omega = 2 * np.pi * f  # angular frequencies for each neuron
    beta = 5  # noise level

    # Filter parameters
    fs = 1/dt
    band_lowcut = f - 2
    band_highcut = f + 2
    lowpass_cutoff = 0.2

    # Solve the ODEs
    x_values, y_values = solve_ode_network(num_steps, dt, a, omega, beta, C, G)

    # Preallocate arrays for filtered signals
    filtered_x = np.zeros_like(x_values)
    # filtered_y = np.zeros_like(y_values)

    # Band-pass filter for each region
    for i in range(78):
        filtered_x[:, i] = bandpass_filter(x_values[:, i], band_lowcut, band_highcut, fs)
        # filtered_y[:, i] = bandpass_filter(y_values[:, i], band_lowcut, band_highcut, fs)
        
    # Calculate the Hilbert transform to get the instantaneous amplitude
    analytical_x = hilbert(filtered_x)
    # analytical_y = hilbert(filtered_y)

    amplitude_env_x = np.abs(analytical_x)
    # amplitude_env_y = np.abs(analytical_y)

    # Smoothening the envelope
    smoothened_env_x = smoothening_envelope(amplitude_env_x, filtered_x)
    # Low-pass filter the amplitude envelope
    analytical_x = hilbert(smoothened_env_x)
    ultra_slow_x = np.zeros_like(smoothened_env_x)
    # ultra_slow_y = np.zeros_like(amplitude_env_y)

    for i in range(78):
        # ultra_slow_x[:, i] = lowpass_filter(amplitude_env_x[:, i], lowpass_cutoff, fs)
        ultra_slow_x[:, i] = lowpass_filter(smoothened_env_x[:, i], lowpass_cutoff, fs)
        # ultra_slow_y[:, i] = lowpass_filter(amplitude_env_y[:, i], lowpass_cutoff, fs)

    return x_values, filtered_x, analytical_x, amplitude_env_x, smoothened_env_x, ultra_slow_x

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