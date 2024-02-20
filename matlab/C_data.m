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

% -------------------------------------------------------------------------

    % creation and use of Euler integration method

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

    beta=55*ones(length(C),1);
    a=0*ones(length(C));
    w=12*2*pi*ones(length(C));
    G=0.5;
    tstep=6000;
    dt=1/100;
    

    


    x=0.5*ones(length(C),1);
    y=0.5*ones(length(C),1);

    [x_val,y_val,time_steps] = Euler2(x,y,tstep);

    

    % ---------------------------------------------------------------------

    % filtering

    

    % bandpass filtering
    function x_val_bandpassed = bandpassing(min,max,x)

        [m,n]=size(x);
    
        x_val_bandpassed = zeros(m,n);

        for i =1:m
            filtered=bandpass(x(i,:),[min max],1000);
            x_val_bandpassed(i,:)=filtered;
    
        end

    end
    

    function x_val_lowpassed = lowpassing(t,x)

        [m,n]=size(x);
    
        x_val_lowpassed = zeros(m,n);
        
        for i =1:m

            filtered=lowpass(x_val(i,:),t,1000);
            x_val_lowpassed(i,:)=filtered;
    
        end

    end

    

    function x_val_envelopes = enveloping(x)
       
        [m,n] = size(x);

        x_val_envelopes=zeros(m,n);

        for i=1:m

            [upper,lower] = envelope(x(i,:));
            x_val_envelopes(i,:) = upper;

        end

    end

    x_val_bandpassed = bandpassing(10,14,x_val);

    
    x_val_lowpassed = lowpassing(0.02,x_val);

    figure(2)
    clf

    subplot(3,1,1)

    plot(time_steps,x_val(45,:))
    hold on
    [upper,lower]=envelope(x_val(45,:));
    plot(time_steps,upper,'linewidth',0.5)
    ylabel('x')

    legend('x(t)','Envelope')

    f_str=num2str(w(45,1)/(2*pi));
    title('bandpass filtering for 10-14Hz with frequency = ', f_str)

    subplot(3,1,2)

    plot(time_steps,x_val_bandpassed(45,:))
    hold on
    [upper_bp,lower]=envelope(x_val_bandpassed(45,:));
    plot(time_steps,upper_bp,'linewidth',0.5)
    ylabel('x')
    legend('x(t)','Envelope')
    

    subplot(3,1,3)
    
    plot(time_steps,x_val_lowpassed(45,:))
    hold on
    [upper_lp,lower]=envelope(x_val_lowpassed(45,:));
    plot(time_steps,upper_lp,'linewidth',0.5)
    ylabel('x')
    xlabel('time steps')
    legend('x(t)','Envelope')

% -------------------------------------------------------------------------

    % hilbert transform (using bandpassed x values)

    % envelope correlation

    bandpassed_envelopes = enveloping(x_val_bandpassed);


    % hilbert transform
    function x_vals_hilbert = Hilbert(x)

        [m,n] = size(x);
        x_vals_hilbert = zeros(m,n);

        for i=1:m

        x_val_hilbert=hilbert(x(i,:));
        inst_phase=angle(x_val_hilbert);

        x_vals_hilbert(i,:) = inst_phase;

        end

    end



    figure(3)
    clf


    x_val_bandpassed = bandpassing(2,6,x_val);
    bandpassed_envelopes = enveloping(x_val_bandpassed);
    envelopes_hilbert = Hilbert(bandpassed_envelopes);
    mean_hilbert = mean(exp(envelopes_hilbert*1i));
    plot(time_steps,mean_hilbert)

    hold on

    x_val_bandpassed = bandpassing(10,14,x_val);
    bandpassed_envelopes = enveloping(x_val_bandpassed);
    envelopes_hilbert = Hilbert(bandpassed_envelopes);
    mean_hilbert = mean(exp(envelopes_hilbert*1i));

    plot(time_steps,mean_hilbert)


    x_val_bandpassed = bandpassing(26,30,x_val);
    bandpassed_envelopes = enveloping(x_val_bandpassed);
    envelopes_hilbert = Hilbert(bandpassed_envelopes);
    mean_hilbert = mean(exp(envelopes_hilbert*1i));

    plot(time_steps,mean_hilbert)

    [M,N]=size(envelopes_hilbert);
    xlabel('Time steps')
    ylabel('Synchronicity')
    title('Change in synhronicity with changing bandpassing range')
    legend('2-6Hz','10-14Hz','26-30Hz')
    xlabel('Time steps')
    ylabel('Synchronicity')
    title('Synchronicity Between brain Regions')

    difference_matrix=zeros(M,M);
    x_val_bandpassed = bandpassing(26,30,x_val);
    bandpassed_envelopes = enveloping(x_val_bandpassed);
    envelopes_hilbert = Hilbert(bandpassed_envelopes);
    mean_hilbert = mean(exp(envelopes_hilbert*1i));

            for M1=1:M
                for M2=1:M

                    difference = mean(cos(envelopes_hilbert(M1)-envelopes_hilbert(M2)));

                    difference_matrix(M1,M2) = difference;

                end
            end


        figure(4)

        imagesc(difference_matrix);
        colorbar;
        xlabel('Brain Region')
        ylabel('Brain Region')
        title('Functional Connectivity Between Different Regions')

    


        








    
    

  
    


end
