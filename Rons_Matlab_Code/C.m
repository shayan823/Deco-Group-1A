function C()

C78 = C_data;
num_steps = 6000;
dt = 1/1000;
P.a = 0;
P.omega = 12*2*pi;
P.beta = 0.5;
P.G = 0.5;
hilbert_envelope_array = [];
angle_array = [];
exp_hilbert_array = [];

% Initialize the arrays to store the x, y values for all neurons
x_values = zeros(num_steps, length(C78));
y_values = zeros(num_steps, length(C78));

% Initialize x and y vectors with zeros
x = ones(1,length(C78))*0.5;
y = ones(1,length(C78))*0.5;

for step = 1:num_steps;
    % Generate a single random noise term for each neuron
    P.noise = randn(1,length(C78))*sqrt(dt);

    % Calculate the dxdt and dydt using the equations provided
    dxdt = (P.a - x.^2 - y.^2).*x - P.omega.*y + (P.G*(diag(C78.*(x - x'))))';
    dydt = (P.a - x.^2 - y.^2).*y + P.omega.*x + (P.G*(diag(C78.*(y - y'))))';

    % Update the x and y values
    x = x+dxdt.*dt+P.beta.*P.noise;
    y = y+dydt.*dt+P.beta.*P.noise;

    % Store the values
    x_values(step,:) = x;
    y_values(step,:) = y;

end;

tstep = 1:6000;

figure(1);
clf;
subplot(2,1,1);
plot(tstep,x_values(:,25));
xlabel("timesteps")
ylabel("x(t)")
hold on;
    
% [yupper,ylower] = envelope(x_values(:,k));
% plot(tstep, yupper);

%% Range Filtering [lower upper]

% Filtering the Original Signal through a 10-14Hz Band
filtered_x_val = bandpass(x_values(:,25), [10 14], 1000);

% Plotting the Filtered Signal (10-14Hz)
subplot(2,1,2);
plot(tstep, filtered_x_val);
xlabel("timesteps")
ylabel("x(t)")
grid on;
hold on;

% Plotting the Upper Envelope on the Subplot (4,1,2)
[yupper,ylower] = envelope(filtered_x_val);
plot(tstep, yupper);

hilbert_envelope = hilbert(yupper);

hilbert_envelope_array = [hilbert_envelope_array hilbert_envelope];

%% Below is Tobias' Method for Calculating CCD and Envelope FC :)
% Bandpassing Envelope Followed By Lowpassing (0.2)
for k = 1:1:length(C78)
    
    filtered_x_val = bandpass(x_values(:,k),[8 12],1000);

    [yupper,ylower] = envelope(filtered_x_val);
    hilbert_envelope = hilbert(yupper);

    hilbert_envelope_array(:,k)= hilbert_envelope;

    low_pass_envelope = lowpass(hilbert_envelope,0.2);
    low_pass_filtered_envelope(:,k) = low_pass_envelope;

end

angle_array = [angle(hilbert_envelope_array)];
exp_phase = exp(1i*angle_array);
phase_sum = abs(sum(exp_phase,2));
R = phase_sum/length(C78);

figure(2);
clf;
plot(tstep,R,'k-');

V=[];
CCD=[];
[m,n] = size(angle_array);

% Envelope FC
figure(3)
clf
envelope_FC=corrcoef(abs(low_pass_filtered_envelope));
envelope_FC=envelope_FC-(diag(envelope_FC).*eye(length(envelope_FC)));
imagesc(envelope_FC)
colormap("turbo");
colorbar;

for time=1:1:m

    for i=1:1:n

        diffs = angle_array(time,i) - angle_array(time,:);
        abs_diff = cos(sqrt(sum(diffs.^2)));
        V(time,i) = abs_diff;

    end
end

for t1 = 1:length(V)
    for t2 = 1:length(V)
        CCD(t1,t2)=(V(t1,:)*V(t2,:)')/(sqrt(sum((V(t1,:).^2)))*sqrt(sum((V(t2,:).^2))'));
    end
end

% CCD
figure(4)
clf
imagesc(CCD);
colormap("turbo");
colorbar;
caxis([0 1]);

end

%% C78
function C78 = C_data;
data = load("C78.mat");
C78 = data.C;
end

