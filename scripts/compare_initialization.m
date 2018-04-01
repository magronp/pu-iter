% Compare different initialization of the iterative procedure
clc; clear all; close all;
test_or_dev = 'Test';
set_settings_puiter;

% Initialization parameters
Nrandom_init = 10;
init = {'Mix','Random','PU'}; Ninit = length(init);
tplot = 255; fspec = 10;

% Define errors and scores
err_rd = zeros(iter_puiter+1,Nrandom_init,Nsongs);
err_pu = zeros(iter_puiter+1,Nsongs);
err_mix = zeros(1,Nsongs);
SDR = zeros(Nsongs,Ninit); SIR = zeros(Nsongs,Ninit); SAR = zeros(Nsongs,Ninit);

for it =1:Nsongs
    
    % Load data
    clc; fprintf('Data %d / %d \n',it,Nsongs);
    num_piece=datavec(it);
    [sm,x,Sm,X,ts,freq] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
    [F,T,J] = size(Sm);
    xe = zeros(J,length(x),Ninit);
    
    % Onset detection
    UN = detect_onset_frames(abs(Sm),Fs,hann(Nw),hop);

    % iterative with random initialization
    for nn=1:Nrandom_init
        [aux,er] = phase_unwrap_ssep(X,Sm,UN,hop,iter_puiter,tplot,0);
        err_rd(:,nn,it)=er(fspec,:);
    end
    Xr=aux; % select one randomly initialized estimate
    
    % Iterative with unwrapping initialization
    [Xu,er]=phase_unwrap_ssep(X,Sm,UN,hop,iter_puiter,tplot,1);
    err_pu(:,it) = er(fspec,:);
 
    % Using the mixture phase
    Xm = abs(Sm).*exp(1i * repmat(angle(X),[1 1 J]));
    err_mix(it) = abs(sum(Xm(fspec,tplot),3)-X(fspec,tplot)).^2;
    
    % Synthesis
    xe(:,:,1)=real(iSTFT(Xm,Nfft,hop,Nw,wtype));
    xe(:,:,2)=real(iSTFT(Xr,Nfft,hop,Nw,wtype));
    xe(:,:,3)=real(iSTFT(Xu,Nfft,hop,Nw,wtype));

    % Score
    for in=1:Ninit
        [sdr,sir,sar] = GetSDR(squeeze(xe(:,:,in)),sm);
        SDR(it,in) = mean(sdr); SIR(it,in) = mean(sir); SAR(it,in) = mean(sar);
    end
   
end

% Save score
save(strcat(metrics_path,'compare_initialization.mat'),'SDR','SIR','SAR','err_pu','err_rd','err_mix');

% % Plot average errors
er_rd = squeeze(mean(err_rd,3));
er_pu = mean(err_pu,2);
er_mi = mean(err_mix);

figure;
semilogy(0:iter_puiter,er_mi*ones(1,iter_puiter+1),'k.',0:iter_puiter,mean(er_rd,2),'b--',0:iter_puiter,er_pu,'r');
xlabel('Iterations','FontSize',16); ylabel('$\mathcal{C}$','interpreter','latex','FontSize',16);