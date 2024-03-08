function mean_correlation

C78 = C_data;
dt = 0.001;
tmax = 300;
sample_rate = 2/dt;
tspan = 0:dt:tmax;
num_steps = length(tspan);
P.a = 0;
P.omega = 12*2*pi;
P.beta = 0.5;
P.G = 0.5;

%Rescale weights as in paper
% C78 = C78*0.2/max(max(C78));

% - Graph Laplacian
Lap = C78 - diag(sum(C78,1));
 
%Plot the structural connectivity
figure(1)
clf
imagesc(C78);
colormap("turbo")
colorbar;
caxis([0 1])
title("Structural Connectivity")

% Initialize the arrays to store the x, y values for all neurons
x_values = zeros(length(C78),num_steps);
y_values = zeros(length(C78),num_steps);

% Initialize x and y vectors
x = ones(length(C78),1)*0.5;
y = ones(length(C78),1)*0.5;

%Bandpass at preferred frequency - for a given node
% node = 10;
% bandpass_signal = bandpass(x_values(node,:), [10 14], sample_rate);
% [amplitude_upper,amplitude_lower] = envelope(bandpass_signal);
% 
% % amplitude_upper = abs(hilbert(real(bandpass_signal)));
% slow_envelope = lowpass(amplitude_upper, 1, sample_rate);

% figure(2);
% clf;
% subplot(3,1,1);
% plot(tspan,x_values(10,:));
% xlabel("timesteps")
% ylabel("x(t)")
% hold on;
% 
% subplot(3,1,2);
% plot(tspan,bandpass_signal);
% xlabel("timesteps")
% ylabel("x(t)")
% hold on;
% plot(tspan,amplitude_upper,'r-');
% 
% subplot(3,1,3);
% plot(tspan,slow_envelope);
% xlabel("timesteps")
% ylabel("x(t)")
% hold on;

beta_values = [0.5, 1.5, 3];
freqs = [2 6; 6 10; 10 14; 14 18; 18 22];

figure(3)
clf
figure(4)
clf

for b = 1:length(beta_values)

    beta=beta_values(b);

        %Euler-Maryuma
    for step = 1:num_steps
    
        % Generate a single random noise term for each neuron
        P.noise = randn(length(C78),1)*sqrt(dt);
    
        % Calculate the dxdt and dydt
        dxdt = (P.a - x.^2 - y.^2).*x - P.omega*y + P.G*Lap*x;
        dydt = (P.a - x.^2 - y.^2).*y + P.omega*x + P.G*Lap*y;
    
        % Update the x and y values
        x = x+dxdt*dt+beta.*P.noise;
        y = y+dydt*dt+beta.*P.noise;
    
        % Store the values
        x_values(:,step) = x;
        y_values(:,step) = y;
    
    end


    M_values = zeros(1, length(freqs));

    for freq=1:length(freqs)

    % now all nodes
    bandpass_signal = bandpass(x_values', freqs(freq,:), sample_rate);
    [amplitude_upper,amplitude_lower] = envelope(bandpass_signal);
    slow_envelope = lowpass(amplitude_upper, 1, sample_rate);
    
    
    hilbert_envelope = hilbert(amplitude_upper);
    angle_envelope = angle(hilbert_envelope);
    
    exp_phase=exp(1i*slow_envelope);
    phase_sum=sqrt(sum(exp_phase.^2,2));
    R=phase_sum/length(C78);
    
    FC = corrcoef(slow_envelope);
    
    %Plot the functional connectivity
    figure(3)
    subplot(length(beta_values),length(freqs),(5*(b-1))+freq)
    imagesc(FC)
    colormap("turbo")
    colorbar;
    % shading interp;
    clim([0 1]);
    title("Functional Connectivity")
    
    % 
    % figure(4)
    % clf
    % correlation_matrix=corrcoef(slow_envelope);
    % correlation_matrix=correlation_matrix-(diag(correlation_matrix).*eye(length(correlation_matrix)));
    % imagesc(abs(correlation_matrix));
    % colorbar;
    % 
    % xlabel('Region')
    % ylabel('Region')
    % title('Envelope FCs')
    
    % V=[];
    % CCD=[];
    % [m,n] = size(slow_envelope);
    % 
    % for time=1:1:m
    % 
    %     for i=1:1:n
    % 
    %         diffs = angle_envelope(time,i)-angle_envelope(time,:);
    %         abs_diff=cos(sqrt(sum(diffs.^2)));
    %         V(time,i)=abs_diff;
    % 
    %     end
    % end
    % 
    % for t1=1:m
    %     for t2=1:m
    %         CCD(t1,t2)=(dot(V(t1,:),V(t2,:)))/(sqrt(sum((V(t1,:).^2)))*sqrt(sum((V(t2,:).^2))));
    %     end
    % end
    % 
    % figure(5)
    % clf
    % imagesc(CCD);
    % colorbar;
    % caxis([0 1])
    % xlabel('Time Step, t_1')
    % ylabel('Time Step, t_2')
    % title('CCD of Brain Regions')

    comparison = corr(C78',FC');
    % figure(6)
    % clf
    % imagesc(comparison)
    % colormap('turbo')
    % colorbar;
    
    M = mean(comparison,'all');
    
    M_values(freq) = M;
    
    freq_ranges = [freqs(:,1), freqs(:,2)];
    
    end

figure(4)
subplot(1,length(beta_values),b)
plot(1:length(M_values), M_values, '-o', 'LineWidth', 2, MarkerFaceColor='auto', DisplayName= ['\beta = ', num2str(P.beta)]);
hold on

end

ylim([0 1])
xlabel('Bandpassing Frequency Range');
ylabel('Mean Correlation');
xticks(1:length(freqs));
xticklabels(cellstr(num2str(freq_ranges)));
xtickangle(45);
grid on;

end

%% C78
function C78 = C_data
data = load("C78.mat");
C78 = data.C;
end
