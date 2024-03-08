function multi_freq_model

function C = C_data;
    data = load("C78.mat");
    C = data.C;
end


C = C_data;
C = C*0.2/max(max(C));
tmax=10;
dt = 0.002;
time_steps=[0:dt:tmax];
num_steps=length(time_steps);
sample_rate=4/dt;

Omegas=[8 12 28]*2*pi;


a = 0;
beta = 0.5;
G = 0.5;



% Initialize the arrays to store the x, y values for all neurons
x_values = zeros(length(C),num_steps);
y_values = zeros(length(C),num_steps);

% Initialize x and y vectors with zeros
x = ones(length(C),1)*0.5;
y = ones(length(C),1)*0.5;

figure(2)
clf
figure(3)
clf
figure(4)
clf

for frequency = 1:length(Omegas)

    omega=Omegas(frequency);
    omega_str=num2str(omega/(2*pi));

    for step = 1:num_steps
    
        % Generate a single random noise term for each neuron
        noise = randn(length(C),1)*sqrt(dt);
    
        % Calculate the dxdt and dydt using the equations provided
        % dxdt = (P.a - x.^2 - y.^2).*x - P.omega*y + (P.G*(diag(C.*(x - x'))))';
        % dydt = (P.a - x.^2 - y.^2).*y + P.omega*x + (P.G*(diag(C.*(y - y'))))';
    
        laplac= C - diag(sum(C,1));
        
        dxdt = (a - x.^2 - y.^2).*x - omega.*y + G*laplac*x;
        dydt = (a - x.^2 - y.^2).*y + omega.*x + G*laplac*y;
    
        % Update the x and y values
        x = x+dxdt*dt+beta*noise;
        y = y+dydt*dt+beta*noise;
    
        % Store the values
        x_values(:,step) = x;
        y_values(:,step) = y;
    
    end
    
    filtered_x_val=bandpass(x_values',[omega-2 omega+2],2000);

    [yupper,ylower]=envelope(filtered_x_val);


    % low_passed_envelope=lowpass(yupper,0.2,2000);
    low_passed_envelope = yupper;
  

    hilbert_envelope = hilbert(yupper);
    angle_envelope = angle(hilbert_envelope);

    instantaneous_phases=angle(low_passed_envelope)
    R=mean(exp(1i*instantaneous_phases),2);
    
    exp_phase=exp(1i*angle_envelope);
    phase_sum=abs(sum(exp_phase,2));
    R=phase_sum/length(C);

    figure(2)

    subplot(2,4,frequency)

    plot(time_steps,R)
    xlabel('Time')
    ylabel('Synchronicity, R(t)')
    title('Synchronicity over time, frequency: ',omega_str )

    figure(3)

    hold on

    subplot(2,4,frequency)

    correlation_matrix=corrcoef(low_passed_envelope);
    correlation_matrix=correlation_matrix-(diag(correlation_matrix).*eye(length(correlation_matrix)));
    imagesc(correlation_matrix);
    colormap('turbo');
    colorbar;
    caxis([-0.1 0.5]);

    xlabel('Region')
    ylabel('Region')
    title('Envelope FCs, frequency: ', omega_str)

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

     subplot(2,4,frequency)

     imagesc(CCD);
     colormap('turbo');
     colorbar;
     caxis([0 1])


    xlabel('Time Step, t_1')
    ylabel('Time Step, t_2')

    title('CCD of Brain Regions, bandpassed: ', omega_str)
end

end
