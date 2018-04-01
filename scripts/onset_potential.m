% Compare different initialization of the iterative procedure
clc; clear all; close all;
test_or_dev = 'Test';
set_settings_puiter;

errp =0:0.1:1; Nerrp = length(errp);

SDRw = zeros(Nsongs,1); SIRw = zeros(Nsongs,1); SARw = zeros(Nsongs,1);
SDRe = zeros(Nsongs,Nerrp); SIRe = zeros(Nsongs,Nerrp); SARe = zeros(Nsongs,Nerrp);
mw = zeros(1,Nsongs);

for it =1:Nsongs
   
    % Source generation
    clc; fprintf('Data %d / %d \n',it,Nsongs);
    num_piece = datavec(it);
    [sm,x,Sm,X,ts,freq] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
    [F,T,J] = size(Sm);
    V = abs(Sm);
    xew = zeros(J,length(x));
    
    % Onset frames detection ans Onset phase (O- and W-)
    UN = detect_onset_frames(abs(Sm),Fs,hann(Nw),hop);
    
    % Iter with PU - mixture onset phase
    Xi_mix = V .* repmat(exp(1i * angle(X)),[1 1 J]);
    Xew = phase_unwrap_ssep(X,Xi_mix,UN,hop,0,1,iter_puiter);
    for j=1:J
       xew(j,:)=real(iSTFT(Xew(:,:,j),Nfft,hop,Nw,wtype));
    end     
    [sdr,sir,sar] = GetSDR(xew,sm);
    SDRw(it) = mean(sdr); SIRw(it) = mean(sir); SARw(it) = mean(sar);
    
    % Adding error to the Oracle
    for ne =1:Nerrp
        fprintf('Error %d / %d \n',ne,Nerrp);
        er = errp(ne);
        Xi_err = V .* exp(1i * ( angle(Sm)+2*pi*er*(randn(F,T,J))));
        % Iter with PU
        Xe_err = phase_unwrap_ssep(X,Xi_err,UN,hop,0,1,iter_puiter);       
        % Synthesis
        xer = zeros(J,length(x));
        for j=1:J
            xer(j,:)=real(iSTFT(Xe_err(:,:,j),Nfft,hop,Nw,wtype));
        end
        %Score
        [sdr,sir,sar] = GetSDR(xer,sm);
        SDRe(it,ne) = mean(sdr); SIRe(it,ne) = mean(sir); SARe(it,ne) = mean(sar);
    end
    
    
    % Compute the error mean and variance (between oracle and Wiener)
    vo = []; vm=[];
    for j=1:J
        vo = [vo angle(Sm(:,UN(j,:)==1,j))];
        vm = [vm angle(X(:,UN(j,:)==1))];
    end
    errwien = abs((vo-vm)./(pi)); errwien = errwien(:);
    mw(it) = mean(errwien);
    
end

% Save results
save(strcat(metrics_path,'onset_potential.mat'),'SDRe','SIRe','SARe','SDRw','SIRw','SARw','errp','mw');

% Plot results
Nerrp = length(errp);
figure;
subplot(1,3,1); plot(errp,mean(SDRw)*ones(1,Nerrp),'k--',errp,mean(SDRe,1),'b'); title('SDR (dB)','fontsize',16); xlabel('\epsilon','fontsize',16);
subplot(1,3,2); plot(errp,mean(SIRw)*ones(1,Nerrp),'k--',errp,mean(SIRe,1),'b'); title('SIR (dB)','fontsize',16); xlabel('\epsilon','fontsize',16);
subplot(1,3,3); plot(errp,mean(SARw)*ones(1,Nerrp),'k--',errp,mean(SARe,1),'b'); title('SAR (dB)','fontsize',16); xlabel('\epsilon','fontsize',16);
hl = legend('Mixture','Oracle+error'); set(hl,'fontsize',14);
