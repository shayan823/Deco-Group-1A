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
C78 = C_data;

%% Plot Nullclines + Trajectories

xmesh = linspace(-1,1,5000);
ymesh = linspace(-1,1,5000);
[x, y] = meshgrid(xmesh,ymesh);

figure(1)
clf
subplot(2,1,1)
hold on

% Plotting Nullclines
contour(x,y,xdot(x,y,P),[0 0],'k','linewidth',2);
contour(x,y,ydot(x,y,P),[0 0],'b','linewidth',2);

% Plotting Trajectory Starting at (0.05, 0.05)
%[t,z] = ode23(@RHS, [t0:0.001:tmax], [0.05 0.05], [], P);

% x = z(:,1);
% y = z(:,2);
% 
% plot(x,y,'g-','LineWidth',1);
% hold on

% Plotting Second Trajectory Starting at Initial Values (at the start of the code)
[t,z] = ode23(@RHS, [t0:0.001:tmax], z0, [], P);

x = z(:,1);
y = z(:,2);

plot(x,y,'r-','LineWidth',1);
xlabel('x')
ylabel('y')
title(['x0 = ', num2str(x0),', y0 = ', num2str(y0),', a = ',num2str(P.a)]);
legend('x-nullcline','y-nullcline','trajectory')
%legend('x-nullcline','y-nullcline','z0 = [0.05 0.05]','z0 = [0.5 0.5]');
hold on

figure(2)
clf
subplot(2,1,1)
plot(t,x,'k','linewidth',2);
xlim([0 20]);
ylim([-2 2]);
xlabel('t(ms)');
ylabel('x');
title(['x0 = ', num2str(x0),', y0 = ', num2str(y0),', a = ',num2str(P.a)]);
hold on

% Creating the second subplot
P.a = -1;

xmesh = linspace(-1,1,700);
ymesh = linspace(-1,1,700);
[x, y] = meshgrid(xmesh,ymesh);

figure(1)
subplot(2,1,2)
hold on

% Plotting Nullclines
contour(x,y,xdot(x,y,P),[0 0],'k','linewidth',2);
contour(x,y,ydot(x,y,P),[0 0],'b','linewidth',2);

% % Plotting Trajectory Starting at (0.05, 0.05)
% [t,z] = ode23(@RHS, [t0:0.001:tmax], [0.05 0.05], [], P);
% 
% x = z(:,1);
% y = z(:,2);
% 
% plot(x,y,'g-','LineWidth',1);
% hold on

% Plotting Second Trajectory Starting at Initial Values (at the start of the code)
[t,z] = ode23(@RHS, [t0:0.001:tmax], z0, [], P);

x = z(:,1);
y = z(:,2);

plot(x,y,'r-','LineWidth',1);
xlabel('x')
ylabel('y')
title(['x0 = ', num2str(x0),', y0 = ', num2str(y0),', a = ',num2str(P.a)]);
legend('x-nullcline','y-nullcline','trajectory')
%legend('x-nullcline','y-nullcline','z0 = [0.05 0.05]','z0 = [0.5 0.5]');
hold on


% %% Plotting x and y Over Time
% 
figure(2)
subplot(2,1,2)
plot(t,x,'k','linewidth',2);
xlim([0 20]);
ylim([-2 2]);
xlabel('t(ms)');
ylabel('x');
title(['x0 = ', num2str(x0),', y0 = ', num2str(y0),', a = ',num2str(P.a)]);
hold on
% 
% figure(3)
% clf
% plot(t,y,'k-','linewidth',1);
% xlim([0 20]);
% ylim([-2 2]);
% xlabel('t(ms)');
% ylabel('y');
% title(['x0 = ', num2str(x0),' y0 = ', num2str(y0),' a = ',num2str(P.a)]);
% hold on
% 
% %% Bifurcation Diagram
% 
% % Bifurcation For -1 < a < 1
% a_range = linspace(-1,1,700);
% 
% x_array = [];
% y_array = [];
% a_array = [];
% 
% for a = a_range
%     P.a = a;
%     [t,z] = ode23(@RHS, [0 20], [0 1], [], P);
% 
%     stored_points = 50;
%     finish = z(end-stored_points:end, :);
% 
%     x_array = [x_array; finish(:, 1)];
%     y_array = [y_array; finish(:, 2)];
%     a_array = [a_array; repmat(a, stored_points+1, 1)];
% 
% end
% 
% figure(4)
% clf
% plot(a_array, x_array,'.', MarkerSize=0.5);
% title('Bifurcation Diagram');
% hold on
% 
% %% Plotting the Structural Connectivity
% 
% figure(5)
% clf
% imagesc(C78);
% colorbar;
% title('Structural Connectivity between Brain Regions');
% xlabel('Brain Area');
% ylabel('Brain Area');
% 
% grid on;
% axis square;
% 
% %% Re-Creating Graphs for Equations 5 and 6
% 
% x = 0.5*ones(length(C78),1);
% y = 0.5*ones(length(C78),1);
% 
% P.beta = 50*ones(length(C78));
% P.a = 0*ones(length(C78));
% P.omega = 2*pi*ones(length(C78));
% P.G = 0.05;
% tstep = 6140;
% dt = 1/100;
% P.noise = normrnd(0,sqrt(dt),[1,length(C78)]);
% 
% [x_val,y_val,time_steps] = Euler(x,y,78,tstep);
% 
% % Plotting the Original Signal
% figure(6)
% clf
% subplot(4,1,1)
% plot(time_steps,x_val,LineWidth=1.5)
% xlim([0 6000])
% ylim([-2 2])
% grid on
% hold on
% 
% % Plotting the Upper Envelope on the Subplot (4,1,1)
% [yupper,ylower] = envelope(x_val);
% plot(time_steps, yupper, LineWidth=2);
% legend('Original Signal','Upper Bound')
% 
% %% Range Filtering [lower upper]
% 
% % Filtering the Original Signal through a 10-14Hz Band
% filtered_x_val = bandpass(x_val, [10 14], 1000);
% 
% % Plotting the Filtered Signal (10-14Hz)
% subplot(4,1,2)
% plot(time_steps, filtered_x_val, LineWidth=1.5);
% xlim([0 6000])
% ylim([-2 2])
% grid on
% hold on
% 
% % Plotting the Upper Envelope on the Subplot (4,1,2)
% [yupper,ylower] = envelope(filtered_x_val);
% plot(time_steps, yupper, LineWidth=2);
% legend('Original Signal Passed Through 10-14Hz Filtering','Upper Envelope')
% 
% %% Lowpass Filtering
% 
% % Lowpass Filtering the Original Signal (<0.2Hz)
% filtered_x_val_2 = lowpass(x_val, 0.2, 1000);
% 
% % Plotting the Filtered Signal (<0.2Hz)
% subplot(4,1,3)
% plot(time_steps, filtered_x_val_2, LineWidth=1.5);
% xlim([0 6000])
% ylim([-2 2])
% grid on
% hold on
% 
% % Plotting the Upper Envelope on the Subplot (4,1,3)
% [yupper,ylower] = envelope(filtered_x_val_2);
% plot(time_steps, yupper, LineWidth=2);
% legend('Lowpass Filtered Original Signal (<0.2Hz)','Upper Envelope')
% 
% %% Highpass Filtering
% 
% % Highpass Filtering the Original Signal (>200Hz)
% filtered_x_val_3 = highpass(x_val, 200, 1000);
% 
% % Plotting the Filtered Signal (>200Hz)
% subplot(4,1,4)
% plot(time_steps, filtered_x_val_3, LineWidth=1.5);
% xlim([0 6000])
% ylim([-2 2])
% grid on
% hold on
% 
% % Plotting the Upper Envelope on the Subplot (4,1,4)
% [yupper,ylower] = envelope(filtered_x_val_3);
% plot(time_steps, yupper, LineWidth=2);
% legend('Highpass Filtered Original Signal (>200Hz)','Upper Envelope')
% 
% %% System with Scalar Input
% 
% x=0.5;
% y=0.5;
% P.beta=0;
% P.omega=2*pi;
% P.G=0;
% P.a=0;
% C78=0;
% 
% tstep=614;
% dt=1/100;
% 
% x_val=[x];
% y_val=[y];
% time_steps=[1];
% 
% for i=2:tstep
% 
%     P.noise = normrnd(0,sqrt(dt));
% 
%     time_steps=[time_steps;i];
% 
%     x = x + dxdt_c(x,y,1)*dt;
%     y = y + dydt_c(x,y,1)*dt;
% 
%     x_val=[x_val;x];
%     y_val=[y_val;y];
% 
% end
% 
% figure(7)
% clf
% plot(time_steps,x_val);
% hold on;
% plot(time_steps,y_val);
% xlim([0 650]);
% 
% %% Bifurcation Diagram with Noise (UNCERTAIN IF THIS IS CORRECT)
% 
% % Creating Arrays for Noise
% x_val_array = [];
% a_array = [];
% 
%    for a = a_range
%     x = 0.5*ones(length(C78),1);
%     y = 0.5*ones(length(C78),1);
%     P.beta = 5*ones(length(C78));
%     P.a = a*ones(length(C78));
%     P.omega = 2*pi*ones(length(C78));
%     P.G = 0.05;
%     tstep = 6140;
%     dt = 1/100;
%     Euler(x,y,1,tstep);
% 
%     stored_points = 50;
%     finish = z(end-stored_points:end, :);
% 
%     x_val_array = [x_val_array; finish(:, 1)];
%     a_array = [a_array; repmat(a, stored_points+1, 1)];
% 
%    end
% 
% % Copying Figure 4 onto Figure 8
% fig4 = figure(4);
% axes4 = findobj(fig4, 'Type', 'axes');
% axesNoTitle = copyobj(axes4, fig4);
% title4 = get(axesNoTitle, 'Title');
% delete(title4);
% fig8 = figure(8);
% copyobj(axes4,fig8);
% hold on
% 
% % Superimposing Noise on Figure 8
% plot(a_array, x_val_array,'.', MarkerSize=0.5,Color='k');
% title('Bifurcation Diagram with Noise');
% hold on

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
%% C78
    function C78 = C_data
    data = load("C78.mat");
    C78 = data.C;
    end
%% Functions for Noise
% dxdt
    function out = dxdt_c(x,y,j)
       out = (P.a(j)-x(j).^2-y(j).^2).*x(j) - P.omega(j).*y(j) + P.G*C78(:,j).*(x-x(j)) + P.beta.*P.noise;
    end
% dydt
    function out = dydt_c(x,y,j)
        out = (P.a(j)-x(j).^2-y(j).^2).*y(j) + P.omega(j).*x(j) + P.G*C78(:,j).*(y-y(j)) + P.beta.*P.noise;
    end
% using entire C78
    function [x_val,y_val,time_steps] = Euler(x,y,j,tstep)
    
        time_steps=[1];

        if j>1
                
                x_val=x(j);
                y_val=y(j);
    
            for i=2:tstep
    
                time_steps=[time_steps;i];
    
                P.noise = normrnd(0,sqrt(dt),[1,length(C78)]);
    
                x = x(j) + dxdt_c(x,y,j)*dt;
                y = y(j) + dydt_c(x,y,j)*dt;
    
                x_val=[x_val;x(j)];
                y_val=[y_val;y(j)];
                
            end
        end
    end
end

 