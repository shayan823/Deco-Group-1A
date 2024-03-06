import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_x_values(x_values, filtered_x, smoothened_env, dt, num_steps):
    # Define the mapping from j values to region names
    region_names = {5: 'Left Frontal Inferior Orbital Area', 17: 'Left Inferior Parietal Area'}
    j = np.array([5, 17])
    time = np.arange(0, num_steps * dt, dt)

    fig, axes = plt.subplots(2, len(j), figsize=(14, len(j)*2))

    for i, region_index in enumerate(j):
        region_name = region_names.get(region_index, f'Region {region_index}')
        axes[0,i].plot(time, x_values[:, region_index], label=f'j={region_index}')
        axes[0,i].set_xlabel('Time (secs)')
        axes[0,i].set_title(f'Original Signal for {region_name}')

        # Filtered signal
        axes[1,i].plot(time, filtered_x[:, region_index], label=f'j={region_index}', color='orange')
        axes[1,i].plot(time, smoothened_env[:, region_index], label=f'j={region_index}', color='green')
        axes[1,i].set_xlabel('Time (secs)')
        axes[1,i].set_title(f'Band-Pass Filtered Signal')

    plt.tight_layout()
    plt.show()