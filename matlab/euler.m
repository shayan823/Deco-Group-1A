function [X,Y] = eulerm(tmax,dt,sample_rate,a,omega,beta,G)

function C = C_data;
    data = load("C78.mat");
    C = data.C;
end


C = C_data;
C = C*0.2/max(max(C));
% tmax=300;
% dt = 0.001;
% sample_rate=1/dt;

time_steps=[0:dt:tmax];
num_steps=length(time_steps);
% 
% 
% P.a = 0;
% P.omega = 24*pi;
% P.beta = 0.5;
% P.G = 0.5;



% Initialize the arrays to store the x, y values for all neurons
x_values = zeros(length(C),num_steps);
y_values = zeros(length(C),num_steps);

% Initialize x and y vectors with zeros
x = ones(length(C),1)*0.5;
y = ones(length(C),1)*0.5;

for step = 1:num_steps

    % Generate a single random noise term for each neuron
    noise = randn(length(C),1)*sqrt(dt);

    % Calculate the dxdt and dydt using the equations provided
    % dxdt = (P.a - x.^2 - y.^2).*x - P.omega*y + (P.G*(diag(C.*(x - x'))))';
    % dydt = (P.a - x.^2 - y.^2).*y + P.omega*x + (P.G*(diag(C.*(y - y'))))';

    laplac= C - diag(sum(C,1));
    
    dxdt = (a - x.^2 - y.^2).*x - omega*y + G*laplac*x;
    dydt = (a - x.^2 - y.^2).*y + omega*x + G*laplac*y;

    % Update the x and y values
    x = x+dxdt*dt+beta*noise;
    y = y+dydt*dt+beta*noise;

    % Store the values
    x_values(:,step) = x;
    y_values(:,step) = y;

end
X=x_values;
Y=y_values;
end
