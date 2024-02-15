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

    function out = dxdt_c(x,y,j,noise)

        % requires vectors x and y, as well as ALL parameters
        % takes jth value of each vector & represents differential equation

        out = (a(j)-x(j).^2-y(j).^2).*x(j) - w(j).*y(j) + G.*sum(C(:,j).*((x-x(j)))) ...
               + beta(j).*noise(j);

    end

    function out = dydt_c(x,y,j,noise)

        % requires vectors x and y as well as ALL parameters to be of same
        % size
        % takes jth value of each vector % represents differential equation
        
        out = (a(j)-x(j).^2-y(j).^2).*y(j) + w(j).*x(j) + G.*sum(C(:,j).*((y-y(j)))) ...
               + beta(j).*noise(j);

    end


    % function [x_val,y_val,time_steps] = Euler(x,y,j,tstep)
    % 
    %     time_steps=[1];
    % 
    %     if j>1
    % 
    %             x_val=x(j);
    %             y_val=y(j);
    % 
    %         for i=2:tstep
    % 
    %             time_steps=[time_steps;i];
    % 
    % 
    %             noise = normrnd(0,sqrt(dt),[1,length(C)]);
    % 
    %             x = x(j) + dxdt_c(x,y,j,noise)*dt;
    %             y = y(j) + dydt_c(x,y,j,noise)*dt;
    % 
    %             x_val=[x_val;x(j)];
    %             y_val=[y_val;y(j)];
    %         end
    %     end
    % end
    % 
    % x=0.5*ones(length(C),1);
    % y=0.5*ones(length(C),1);
    % 
    % beta=20*ones(length(C));
    % a=0*ones(length(C));
    % w=2*pi*ones(length(C));
    % G=0.0;
    % tstep=6000;
    % dt=1/100;
    % noise = normrnd(0,sqrt(dt),[1,length(C)]);
    % 
    % [x_val,y_val,time_steps] = Euler(x,y,30,tstep);
    % 
    %  figure(2)
    %  plot(time_steps,x_val)


% -------------------------------------------------------------------------


    % beta=50*ones(length(C));
    % a=0*ones(length(C));
    % w=2*pi*ones(length(C));
    % G=0.0;
    % tstep=6000;
    % dt=1/100;
    % 
    % noise = normrnd(0,sqrt(dt),[1,length(C)]);
    % 
    % 
    % euler_xvals=[];
    % euler_yvals=[];
    % 
    % 
    % 
    % for i = 2:length(C)
    % 
    % 
    %     x=0.5*ones(length(C),1);
    %     y=0.5*ones(length(C),1);
    % 
    % 
    %      [x_val,y_val,time_steps] = Euler(x,y,i,tstep);
    % 
    % 
    %     euler_xvals=[euler_xvals x_val];
    %     euler_yvals=[euler_yvals y_val];
    % 
    % end

     function [x_val,y_val,time_steps] = Euler2(x,y,tstep)
                
                x_val=x;
                y_val=y;
                time_steps=1;

            for i=2:tstep
                
                time_steps=[time_steps;i];

                noise = normrnd(0,sqrt(dt),[length(C),1]);

                for j=1:length(x)

                   x(j) = x(j) + dxdt_c(x,y,j,noise)*dt;
                   y(j) = y(j) + dydt_c(x,y,j,noise)*dt;
                
                end

                x_val=[x_val x];
                y_val=[y_val y];

            end
     end

    beta=50*ones(length(C),1);
    a=0*ones(length(C));
    w=2*pi*ones(length(C));
    G=0.05;
    tstep=6000;
    dt=1/100;
    


    x=0.5*ones(length(C),1);
    y=0.5*ones(length(C),1);

    [x_val,y_val,time_steps] = Euler2(x,y,tstep);

    figure(3)
    clf

    plot(time_steps,x_val(43,:))
    hold on
    % plot(time_steps,x_val(47,:))
    

  
    


end
