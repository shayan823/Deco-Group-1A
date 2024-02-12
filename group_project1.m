function group_project1()
 
%% Initialise and Declare Parameter Values
P.a = 0;
P.omega = 2*pi;

z = [];
t = [];
x0 = .5;
y0 = .5;
z0 = [x0; y0];
t0 = 0;
tmax = 20;

%% Plot phase planes

xmesh = linspace(-1,1,700);
ymesh = linspace(-1,1,700);
[x, y] = meshgrid(xmesh,ymesh);

figure(1)
clf
hold on

contour(x,y,xdot(x,y,P),[0 0],'k','linewidth',2)
contour(x,y,ydot(x,y,P),[0 0],'b','linewidth',2)
legend('x-nullcline','y-nullcline');

[t,z] = ode23(@RHS, [t0:0.001:tmax], z0, [], P)

x = z(:,1);
y = z(:,2);

figure(1)
plot(x,y,'r-','LineWidth',1)
hold on

% zss = fsolve(@(z) SteadyState(z,P), [-0.1;0.1])
% xss = zss(1);
% yss = zss(2);
% figure (1)
% plot(xss,yss,'*g','linewidth',1)

figure(2)
clf
plot(t,x,'k','linewidth',1)
xlim([0 20])
ylim([-2 2])
xlabel('t(ms)')
ylabel('x')
title(['x0 = ', num2str(x0),' y0 = ', num2str(y0),' a = ',num2str(P.a)])

figure(3)
clf
plot(t,y,'k','linewidth',1)
xlim([0 20])
ylim([-2 2])
xlabel('t(ms)')
ylabel('y')
title(['x0 = ', num2str(x0),' y0 = ', num2str(y0),' a = ',num2str(P.a)])

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
    
%% RHS(t,y,P)    
    function RHS = RHS(t,z,P)
        x=z(1);
        y=z(2);
        RHS = [xdot(x,y,P);ydot(x,y,P)];
    end
end
 