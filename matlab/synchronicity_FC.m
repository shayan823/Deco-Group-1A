function synchronicity_FC;

function C = C_data;
    data = load("C78.mat");
    C = data.C;
end


C = C_data;
C = C*0.2/max(max(C));
tmax=3.2;
dt = 1/1000;
time_steps=[0:dt:tmax];
num_steps=length(time_steps);

P.a = 0;
P.omega = 24*pi;
P.beta = 0.5;
P.G = 0.2;



% Initialize the arrays to store the x, y values for all neurons
x_values = zeros(length(C),num_steps);
y_values = zeros(length(C),num_steps);

% Initialize x and y vectors with zeros
x = ones(length(C),1)*0.5;
y = ones(length(C),1)*0.5;

for step = 1:num_steps

    % Generate a single random noise term for each neuron
    P.noise = randn(length(C),1)*sqrt(dt);

    % Calculate the dxdt and dydt using the equations provided
    % dxdt = (P.a - x.^2 - y.^2).*x - P.omega*y + (P.G*(diag(C.*(x - x'))))';
    % dydt = (P.a - x.^2 - y.^2).*y + P.omega*x + (P.G*(diag(C.*(y - y'))))';

    laplac= C - diag(sum(C,1));
    
    dxdt = (P.a - x.^2 - y.^2).*x - P.omega*y + P.G*laplac*x;
    dydt = (P.a - x.^2 - y.^2).*y + P.omega*x + P.G*laplac*y;

    % Update the x and y values
    x = x+dxdt*dt+P.beta*P.noise;
    y = y+dydt*dt+P.beta*P.noise;

    % Store the values
    x_values(:,step) = x;
    y_values(:,step) = y;

end

freqs = [2 6;10 14; 26 30];

figure(3)
clf

figure(2)
clf

 figure(4)
 clf

for freq=1:length(freqs)
    
    frequency=num2str(freqs(freq,:));

    filtered_x_val=bandpass(x_values',freqs(freq,:),1000);

    [yupper,ylower]=envelope(filtered_x_val);


    low_passed_envelope=lowpass(yupper,0.2);

    hilbert_envelope = hilbert(low_passed_envelope);
    angle_envelope = angle(hilbert_envelope);

    % instantaneous_phases=angle(low_passed_envelope)
    % R=mean(exp(1i*instantaneous_phases),2);
    
    exp_phase=exp(1i*low_passed_envelope);
    phase_sum=sqrt(sum(exp_phase.^2,2));
    R=phase_sum/length(C);
    
    figure(2)

    subplot(1,3,freq)
    
    plot(time_steps,R)
    xlabel('Time')
    ylabel('Synchronicity, R(t)')
    title('Synchronicity over time, bandpassed: ',frequency )

    figure(3)
    
    hold on

    subplot(1,3,freq)

    correlation_matrix=corrcoef(low_passed_envelope);
    correlation_matrix=correlation_matrix-(diag(correlation_matrix).*eye(length(correlation_matrix)));
    imagesc(abs(correlation_matrix));
    colorbar;

    xlabel('Region')
    ylabel('Region')
    title('Envelope FCs, bandpassed: ', frequency)

    V=[];
    CCD=[];
    [m,n] = size(low_passed_envelope);

    for time=1:1:m
    
        for i=1:1:n
    
            diffs = angle_envelope(time,i)-angle_envelope(time,:);
            abs_diff=cos(sqrt(sum(diffs.^2)));
            V(time,i)=abs_diff;
    
        end
    end

    for t1=1:m
        for t2=1:m
            
            CCD(t1,t2)=(dot(V(t1,:),V(t2,:)))/(sqrt(sum((V(t1,:).^2)))*sqrt(sum((V(t2,:).^2))));
        end
        % time=t1
    end
    
     figure(4)
     
     subplot(1,3,freq)
    
     imagesc(CCD);
     colorbar;
     caxis([0 1])
    
    
    xlabel('Time Step, t_1')
    ylabel('Time Step, t_2')
    
    title('CCD of Brain Regions, bandpassed: ', frequency)


end



end
