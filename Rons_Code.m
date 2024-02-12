function Rons_Code()
 
%% Initialise and Declare Parameter Values

P.a = 0;
P.omega = 2*pi;

z = []; 
t = [];
x0 = 0.5;
y0 = 0.5;
z0 = [x0; y0];
t0 = 0;
tmax = 20;

%% Plot Nullclines + Trajectories

xmesh = linspace(-1,1,700);
ymesh = linspace(-1,1,700);
[x, y] = meshgrid(xmesh,ymesh);

figure(1)
clf
hold on

contour(x,y,xdot(x,y,P),[0 0],'k','linewidth',2);
contour(x,y,ydot(x,y,P),[0 0],'b','linewidth',2);

[t,z] = ode23(@RHS, [t0:0.001:tmax], [0.05 0.05], [], P);

x = z(:,1);
y = z(:,2);

plot(x,y,'g-','LineWidth',1);
hold on

[t,z] = ode23(@RHS, [t0:0.001:tmax], z0, [], P);

x = z(:,1);
y = z(:,2);

plot(x,y,'r-','LineWidth',1);
legend('x-nullcline','y-nullcline','z0 = [0.05 0.05]','z0 = [0.5 0.5]');
hold on

%% Plotting x and y Over Time

figure(2)
clf
plot(t,x,'k-','linewidth',1);
xlim([0 20]);
ylim([-2 2]);
xlabel('t(ms)');
ylabel('x');
title(['x0 = ', num2str(x0),' y0 = ', num2str(y0),' a = ',num2str(P.a)]);
hold on

figure(3)
clf
plot(t,y,'k-','linewidth',1);
xlim([0 20]);
ylim([-2 2]);
xlabel('t(ms)');
ylabel('y');
title(['x0 = ', num2str(x0),' y0 = ', num2str(y0),' a = ',num2str(P.a)]);
hold on

%% Bifurcation Diagram

a_range = linspace(-1,1,700);

x_array = [];
y_array = [];
a_array = [];

for a = a_range
    P.a = a;
    [t,z] = ode23(@RHS, [0 20], [0 1], [], P);

    stored_points = 50;
    finish = z(end-stored_points:end, :);

    x_array = [x_array; finish(:, 1)];
    y_array = [y_array; finish(:, 2)];
    a_array = [a_array; repmat(a, stored_points+1, 1)];

end

figure(4)
clf
plot(a_array, x_array,'.', MarkerSize=0.5)
title('Bifurcation Diagram')

%% Function Definitions
%% Xdot(X,Y,P)    
    function xdot = xdot(x,y,P)
        xdot = (P.a-x.^2-y.^2).*x-(P.omega.*y);
     end

%%  Ydot(X,Y,P)  
    function ydot = ydot(x,y,P)
        ydot = (P.a-x.^2-y.^2).*y+(P.omega.*x);
    end
%% SteadyState(y,P)    
    function SteadyState = SteadyState(y,P)
        x = z(1);
        y = z(2);
        SteadyState = [xdot(x,y,P);ydot(x,y,P)];
    end
    
%% RHS(t,z,P)    
    function RHS = RHS(t,z,P)
        x=z(1);
        y=z(2);
        RHS = [xdot(x,y,P);ydot(x,y,P)];
    end
end
 