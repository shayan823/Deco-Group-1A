function MEG

        
    % %EM Euler-Maruyama method on linear SDE
    % %
    % % SDE is dX = lambda*X dt + mu*X dW, X(0) = Xzero,
    % % where lambda = 2, mu = 1 and Xzero = 1.
    % %
    % % Discretized Brownian path over [0,1] has dt = 2^(-8).
    % % Euler-Maruyama uses timestep R*dt.
    % randn('state',100)
    % lambda = 2; mu = 1; Xzero = 1; % problem parameters
    % T = 1; N = 2^11; dt = 1/N;
    % dW = sqrt(dt)*randn(1,N); % Brownian increments
    % W = cumsum(dW); % discretized Brownian path
    % Xtrue = Xzero*exp((lambda-0.5*mu^2)*([dt:dt:T])+mu*W);
    % 
    % figure(3)
    % 
    % tt=[0:dt:T]
    % plot(tt,[Xzero,Xtrue],'m-'), hold on
    % R = 4; Dt = R*dt; L = N/R; % L EM steps of size Dt = R*dt
    % Xem = zeros(1,L); % preallocate for efficiency
    % Xtemp = Xzero;
    % 
    % for j = 1:L
    % Winc = sum(dW(R*(j-1)+1:R*j));
    % Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winc;
    % Xem(j) = Xtemp;
    % end
    % 
    % Tt=[0:Dt:T];
    % Xems=[Xzero,Xem];
    % plot(Tt,Xems,'r--*'), hold off
    % xlabel('t','FontSize',12)
    % ylabel('X','FontSize',16,'Rotation',0,'HorizontalAlignment','right')
    % emerr = abs(Xem(end)-Xtrue(end));

    % ODE solutions

    function out = functions(t,p)

        %dx/dt and dy/dt

    x = p(1);
    y = p(2);

    out = [ (a-x.^2-y.^2).*x-w.*y;
            (a-x.^2-y.^2).*y+w.*x];
    end


    % trial

    a = 0.5; % HOPF bifurcations at a=0, when a<0 no oscillations, when a>=0 oscillations
    ws = [4:4:24]; % f=w/2pi => w = f*2pi

    figure(1)
    clf

    figure(2)
    clf

    for i=1:length(ws)

        % with changing frequencies

        w=ws(i)*2*pi;

        figure(1)

        subplot(3,2,i)

        [t,x]=ode45(@functions,linspace(0,1.5,1000),[sqrt(a) sqrt(a)]);

        Xx=x(:,1);
        plot(t,Xx,'linewidth',1)
        hold on

        plot(t,x(:,2),'linewidth',1)
        legend('x','y')

        w=num2str(w);
        title('w = ', w)
        xlabel('t')
        ylabel('x(t)/y(t)')

    end

    as=linspace(-3,2,6);
    w=2*pi;

    for j = 1:length(as)

        % with changing a - shows bifurcation

        a=as(j);

        figure(2)

        subplot(3,2,j)

        [t,x]=ode45(@functions,linspace(0,2,1000),[sqrt(a)+1 sqrt(a)+1]);
        plot(t,x(:,1),'linewidth',1)
        hold on

        plot(t,x(:,2),'linewidth',1)
        % legend('x','y')

        a=num2str(a);
        title('a = ', a)
        xlabel('t')
        ylabel('x(t)/y(t)')


    end



    % using euler's method 

    function out = dxdt(x,y)

        out = (a-x^2-y^2)*x - w*y + beta*noise;

    end

    function out = dydt(x,y)

        out = (a-x^2-y^2)*y + w*x + beta*noise;

    end

    % noisy solutions with changing a

    as=[-2,-1,0,1];

    figure(4)
        clf

for j=1:length(as)

    % initial conditions

    w=2*pi;
    a=as(j);
    tsteps=614;
    dt=1/100;
    beta=50;
    noise=[];

    x=0.5;
    y=0.5;

    x_val=[x];
    y_val=[y];

    time_steps=[1];

    for i=2:tsteps

        % 1st time step already considered, move to next
        % continue from there

        time_steps=[time_steps;i];

        % generate randomly generated, normally distributed value

        noise = normrnd(0,sqrt(dt));

       % euler method

        x = x + dxdt(x,y)*dt;
        y = y + dydt(x,y)*dt;

        % add new value to list of x,y values

        x_val=[x_val;x];
        y_val=[y_val;y];

    end

    % plotting

    figure(4)

    subplot(2,2,j)

    plot(time_steps,x_val)

    hold on

    plot(time_steps,y_val)

    a=num2str(a)


    legend('dx/dt','dy/dt')
    xlabel('time steps')
    ylabel('x(t) / y(t)')
    title('ODE solutions when a = ', a)

end

% normal ODE solutions
% same as previous but without noise

figure(5)
    clf

for j=1:length(as)

    % initial conditions

    w=2*pi;
    a=as(j);
    tsteps=614;
    dt=1/100;
    beta=0;
    noise=[];

    x=0.5;
    y=0.5;

    x_val=[x];
    y_val=[y];

    time_steps=[1];

    for i=2:tsteps

        % 1st time step considered, move to next
        % continue from there

        time_steps=[time_steps;i];

        % generate randomly generated, normally distributed variable

        noise = normrnd(0,sqrt(dt));

        % euler method

        x = x + dxdt(x,y)*dt;
        y = y + dydt(x,y)*dt;

        % add updated value to x,y values

        x_val=[x_val;x];
        y_val=[y_val;y];

    end

    % plot

    figure(5)

    subplot(2,2,j)

    plot(time_steps,x_val)

    hold on

    plot(time_steps,y_val)

    a=num2str(a);


    legend('dx/dt','dy/dt')
    xlabel('time steps')
    ylabel('x(t) / y(t)')
    title('ODE solutions when a = ', a)
end





% bifurcation diagram
    
    function out = system_equ(p)

        x = p(1);
        y = p(2);

        out = [(a-x.^2-y.^2).*x - w*y + beta*noise;
               (a-x.^2-y.^2).*y + w*x + beta*noise;];
    end

    
    
    function out = Jacob(x,y)
    
        out = [a-3*x^2-y^2, -2*x*y-w;
               -2*x*y+w, a-x^2-3*y^2];
    end


as=linspace(-3,3,1000);
w=2*pi;
dt=1/100;
beta=0;
eigenvalues_val=[];
eigenvalues_val2=[];

options = optimoptions('fsolve', 'Display', 'off');

for p=1:length(as)

    a=as(p);
    noise = normrnd(0,sqrt(dt));

    system_eq_handle = @(p) system_equ(p);
    ss = fsolve(system_eq_handle, [0.5, 0.5], options);

    Jacobian_init=Jacob(ss(1),ss(2));
    
    eigenvalues=eig(Jacobian_init);

    eigenvalues_val=[eigenvalues_val,eigenvalues(1)];
    eigenvalues_val2=[eigenvalues_val2,eigenvalues(2)];


end


figure(6)
clf
plot(real(eigenvalues_val),imag(eigenvalues_val))
hold on
plot(real(eigenvalues_val2),imag(eigenvalues_val2))

   
    
end
