from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
        dxdt = (a - x**2 - y**2)*x - omega*y + G * np.dot(C, (x - x[:, None]).T).diagonal() + beta * noise
        dydt = (a - x**2 - y**2)*y + omega*x + G * np.dot(C, (y - y[:, None]).T).diagonal() + beta * noise
        
        # Update the x and y values
        x += dxdt * dt
        y += dydt * dt
        
        # Store the values
        x_values[step, :] = x
        y_values[step, :] = y

    
    return x_values, y_values