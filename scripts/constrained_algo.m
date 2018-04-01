% Compare PU-Iter using with or without adding an extra constraint
clc; clear all; close all;
test_or_dev = 'Test';
set_settings_puiter;

% Constrained algo
diff_constraints = 2;           % either constraint or not (= classical PU-Iter)
Sigma = 10.^(-3:0.5:3);
Ns=length(Sigma);

SDR = zeros(Ns,diff_constraints,Nsongs); SIR = zeros(Ns,diff_constraints,Nsongs); SAR = zeros(Ns,diff_constraints,Nsongs);

for it =1:Nsongs
    
    % Source generation
    num_piece=it;
    [sm,x,Sm,X,ts,freq] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
    [F,T,J] = size(Sm);
    Yrand = abs(Sm) .* exp(1i * 2*pi*rand(size(Sm)));       % test for a random initialization
    xe = zeros(J,length(x),diff_constraints);
    UN = detect_onset_frames(abs(Sm),Fs,hann(Nw),hop);
    
    for ss=1:Ns
        sigma=Sigma(ss);
        clc; fprintf('Data %d / %d \n Sigma %d / %d \n',it,Nsongs,ss,Ns);
        
        % Source separation
        Xr = phase_unwrap_ssep_constrained(X,Sm,UN,hop,iter_puiter,sigma,0,Yrand);
        Xu = phase_unwrap_ssep_constrained(X,Sm,UN,hop,iter_puiter,sigma,1);

        % Synthesis
        for j=1:J
            xe(j,:,1)=real(iSTFT(Xr(:,:,j),Nfft,hop,Nw,wtype));
            xe(j,:,2)=real(iSTFT(Xu(:,:,j),Nfft,hop,Nw,wtype));
        end

        % Score
        for co=1:diff_constraints
            [sdr,sir,sar] = GetSDR(squeeze(xe(:,:,co)),sm);
            SDR(ss,co,it) = mean(sdr); SIR(ss,co,it) = mean(sir); SAR(ss,co,it) = mean(sar);
        end
    end
end

save(strcat(metrics_path,'constrained_algo.mat'),'SDR','SIR','SAR','Sigma');

SDR_av = mean(SDR,3); SIR_av = mean(SIR,3); SAR_av = mean(SAR,3);
figure;
subplot(3,1,1); semilogx(Sigma,SDR_av(:,1),'b+-'); hold on; semilogx(Sigma,SDR_av(:,2),'r*-'); hold on; semilogx(Sigma,ones(1,Ns)*SDR_av(1,2),'k'); hold on; semilogx(Sigma,ones(1,Ns)*SDR_av(end,2),'k--'); aa=axis; axis([Sigma(1) Sigma(end) aa(3) aa(4)]); ha=xlabel('\sigma'); set(ha,'FontSize',16); ha=ylabel('SDR (dB)'); set(ha,'FontSize',16);
subplot(3,1,2); semilogx(Sigma,SIR_av(:,1),'b+-'); hold on; semilogx(Sigma,SIR_av(:,2),'r*-'); hold on; semilogx(Sigma,ones(1,Ns)*SIR_av(1,2),'k'); hold on; semilogx(Sigma,ones(1,Ns)*SIR_av(end,2),'k--'); aa=axis; axis([Sigma(1) Sigma(end) aa(3) aa(4)]); ha=xlabel('\sigma'); set(ha,'FontSize',16); ha=ylabel('SIR (dB)'); set(ha,'FontSize',16);
subplot(3,1,3); semilogx(Sigma,SAR_av(:,1),'b+-'); hold on; semilogx(Sigma,SAR_av(:,2),'r*-'); hold on; semilogx(Sigma,ones(1,Ns)*SAR_av(1,2),'k'); hold on; semilogx(Sigma,ones(1,Ns)*SAR_av(end,2),'k--'); aa=axis; axis([Sigma(1) Sigma(end) aa(3) aa(4)]); ha=xlabel('\sigma'); set(ha,'FontSize',16); ha=ylabel('SAR (dB)'); set(ha,'FontSize',16);
ha=legend('Random','PU'); set(ha,'FontSize',14);
