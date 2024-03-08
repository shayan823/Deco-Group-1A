function CCD(x_values,freq,sample_rate)

    filtered_x_val=bandpass(x_values',freq,sample_rate);
    
        [yupper,ylower]=envelope(filtered_x_val);
    
    
        low_passed_envelope=lowpass(yupper,0.2,sample_rate);
        %  low_passed_envelope = yupper;
      
    
        hilbert_envelope = hilbert(low_passed_envelope);
        angle_envelope = angle(hilbert_envelope);

    V=[];
    CCD=[];
    [m,n] = size(angle_envelope);

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
        % print(t1, 'done')
    end
    
     figure(4)
     
     % subplot(length(P.beta),length(freq),freq)
    
     imagesc(CCD);
     colormap('turbo');
     colorbar;
     caxis([0 1])
    
    
    xlabel('Time Step, t_1')
    ylabel('Time Step, t_2')
    
    title('CCD of Brain Regions, bandpassed: ', frequency)


end
