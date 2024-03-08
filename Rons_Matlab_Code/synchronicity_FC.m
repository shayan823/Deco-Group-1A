function synchronicity_FC

function C78 = C_data
    data = load("C78.mat");
    C78 = data.C;
end

C78 = C_data;
dt = 0.001;
tmax = 10;
sample_rate = 0.5/dt;
tspan = 0:dt:tmax;
num_steps = length(tspan);
P.a = 0;
P.omega = 12*2*pi;
P.beta = 0.5;
P.G = 0.5;

% Initialize the arrays to store the x, y values for all neurons
x_values = zeros(length(C78),num_steps);
y_values = zeros(length(C78),num_steps);

% Initialize x and y vectors with zeros
x = ones(length(C78),1)*0.5;
y = ones(length(C78),1)*0.5;

for step = 1:num_steps

    % Generate a single random noise term for each neuron
    P.noise = randn(length(C78),1)*sqrt(dt);

    laplac= C78 - diag(sum(C78,1));

    dxdt = (P.a - x.^2 - y.^2).*x - P.omega*y + P.G*laplac*x;
    dydt = (P.a - x.^2 - y.^2).*y + P.omega*x + P.G*laplac*y;

    % Update the x and y values
    x = x+dxdt*dt+P.beta*P.noise;
    y = y+dydt*dt+P.beta*P.noise;

    % Store the values
    x_values(:,step) = x;
    y_values(:,step) = y;

end

    filtered_x_val=bandpass(x_values',[10 14],sample_rate);

    [yupper,ylower]=envelope(filtered_x_val);

    hilbert_envelope = hilbert([yupper,ylower]);
    angle_envelope = angle(hilbert_envelope);
    
    figure(1)
    clf
    correlation_matrix=corrcoef([yupper,ylower]);
    correlation_matrix=correlation_matrix-(diag(correlation_matrix).*eye(length(correlation_matrix)));
    imagesc(abs(correlation_matrix));
    colorbar;

    xlabel('Region')
    ylabel('Region')
    title('Envelope FCs')

    V=[];
    CCD=[];
    [m,n] = size([yupper,ylower]);

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
    end

    figure(2)
    imagesc(CCD);
    colorbar;
    caxis([0 1])
    xlabel('Time Step, t_1')
    ylabel('Time Step, t_2')
    title('CCD of Brain Regions')

end