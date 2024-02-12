function C_data

    data=load("C90.mat");
    C=data.C;
    
    figure(1)
    clf
    
    imagesc(C);
    colorbar;
    
    title('Structural Connectivity between Brain Regions')
    xlabel('Brain Area')
    ylabel('Brain Area')
     
    grid on;
    axis square;

% -------------------------------------------------------------------------

    function out = dxdt_c(x,y,j)

        % requires vectors x and y, as well as ALL parameters
        % takes jth value of each vector & represents differential equation

        out = (a(j)-x(j).^2-y(j).^2).*x(j) - w(j).*y(j) + G*C(:,j).*(x-x(j)) ...
               + beta.*noise;

    end

    function out = dydt_c(x,y,j)

        % requires vectors x and y as well as ALL parameters to be of same
        % size
        % takes jth value of each vector % represents differential equation
        
        out = (a(j)-x(j).^2-y(j).^2).*y(j) + w(j).*x(j) + G*C(:,j).*(y-y(j)) ...
               + beta.*noise;

    end


    



    % Using entire C

    function [x_val,y_val,time_steps] = Euler(x,y,j,tstep)

        time_steps=[1];

        if j>1
                
                x_val=x(j);
                y_val=y(j);

            for i=2:tstep

                time_steps=[time_steps;i];

                
                noise = normrnd(0,sqrt(dt),[1,length(C)]);

                x = x(j) + dxdt_c(x,y,j)*dt;
                y = y(j) + dydt_c(x,y,j)*dt;

                x_val=[x_val;x(j)];
                y_val=[y_val;y(j)];
            end
        end
    end

    x=0.5*ones(length(C),1);
    y=0.5*ones(length(C),1);

    beta=5*ones(length(C));
    a=0*ones(length(C));
    w=2*pi*ones(length(C));
    G=0.05;
    tstep=6140;
    dt=1/100;
    noise = normrnd(0,sqrt(dt),[1,length(C)]);

    [x_val,y_val,time_steps] = Euler(x,y,81,tstep);

     figure(3)
     plot(time_steps,x_val)



     % system with scalar input

    x=0.5;
    y=0.5;
    beta=0;
    w=2*pi;
    G=0;
    a=0;
    C=0; % - un-comment to view

    tstep=614;
    dt=1/100;

    x_val=[x];
    y_val=[y];
    time_steps=[1];

    for i=2:tstep


        noise = normrnd(0,sqrt(dt));

        time_steps=[time_steps;i];

        x = x + dxdt_c(x,y,1)*dt;
        y = y + dydt_c(x,y,1)*dt;

        x_val=[x_val;x];
        y_val=[y_val;y];

    end

    x_val
    figure(2)
    clf

    plot(time_steps,x_val)
    hold on
    plot(time_steps,y_val)




end
