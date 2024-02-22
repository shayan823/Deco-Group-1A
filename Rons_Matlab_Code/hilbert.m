function [x_values,y_values] = hilbert()

C78 = C_data;
num_steps = 6000;
dt = 1/100;
P.a = 0;
P.omega = 2*pi;
P.beta = 5;
P.G = 0.5;
hilbert_envelope_array = [];

% Initialize the arrays to store the x, y values for all neurons
x_values = zeros(num_steps, 78);
y_values = zeros(num_steps, 78);

% Initialize x and y vectors with zeros
x = ones(1,78)*0.5;
y = ones(1,78)*0.5;

for step = 1:num_steps
    % Generate a single random noise term for each neuron
    P.noise = randn(1,78)*sqrt(dt);

    % Calculate the dxdt and dydt using the equations provided
    dxdt = (P.a - x.^2 - y.^2).*x - P.omega.*y + (P.G*(diag(C78.*(x - x'))))';
    dydt = (P.a - x.^2 - y.^2).*y - P.omega.*x + (P.G*(diag(C78.*(y - y'))))';

    % Update the x and y values
    x = x+dxdt.*dt+P.beta.*P.noise;
    y = y+dydt.*dt+P.beta.*P.noise;

    % Store the values
    x_values(step,:) = x;
    y_values(step,:) = y;

end

tstep = 1:6000;

figure(1);
clf

for k=1:1:length(C78)

    subplot(2,1,1);
    plot(tstep,x_values(:,k));
    hold on
    
    [yupper,ylower] = envelope(x_values(:,k));
    plot(tstep, yupper);
    
    %% Range Filtering [lower upper]
    
    % Filtering the Original Signal through a 10-14Hz Band
    filtered_x_val = bandpass(x_values(:,k), [10 14], 1000);
    
    % Plotting the Filtered Signal (10-14Hz)
    subplot(2,1,2);
    plot(tstep, filtered_x_val);
    grid on
    hold on
    
    % Plotting the Upper Envelope on the Subplot (4,1,2)
    [yupper,ylower] = envelope(filtered_x_val);
    plot(tstep, yupper);

    hilbert_envelope = hilbert(yupper);

    hilbert_envelope_array = [hilbert_envelope_array hilbert_envelope];
    
    end

end

%% C78
function C78 = C_data
data = load("C78.mat");
C78 = data.C;
end
