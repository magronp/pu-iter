% Compare Griffin Lim algorithm and PU for a blind phase recovery task
clc; clear all; close all;
test_or_dev = 'Test';
set_settings_puiter;

% Window lengths
wlength = 2.^(9:14);
LNw = length(wlength);

iter_GL = 1;
direc = strcat(audio_path,'gl_vs_pu/');
SDR = zeros(LNw,4,Nsongs); icons = zeros(LNw,4,Nsongs);

for it=1:Nsongs
    for nw=1:LNw
        fprintf('Iteration %d / %d \n Wind length %d / %d \n',it,Nsongs,nw,LNw);

        %%% Load data
        num_piece = it;
        Nw = wlength(nw);
        Nfft = Nw; hop = Nw/4;
        [sm,x,Sm,X,ts,freq] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
        [F,T,J] = size(Sm);
        V = abs(X);

        %%% NMF
        W_ini = rand(F,K); H_ini = rand(K,T);
        [W,H] = NMF(V,W_ini,H_ini,iter_nmf,1,0);
        Vapprox = W*H;
        Xapprox = Vapprox.*exp(1i*angle(X));
        
        %%% Phase recovery
        
        % Onset detections
        UN_oracle = detect_onset_frames(V,Fs,hann(Nw),hop); delta_oracle = repmat(UN_oracle,[F 1]); delta_oracle = 1-delta_oracle;
        UN_nmf = detect_onset_frames(Vapprox,Fs,hann(Nw),hop); delta_nmf = repmat(UN_nmf,[F 1]); delta_nmf = 1-delta_nmf;

        % Corrupt with random phase - oracle - magnitude
        phir = 2*pi*rand(F,T); phir(:,UN_oracle==1) = angle(X(:,UN_oracle==1));
        XC_oracle = V .* exp(1i * phir);
        xc_oracle = real(iSTFT(XC_oracle,Nfft,hop,Nw,wtype));

        % Corrupt with random phase - oracle - magnitude
        phir = 2*pi*rand(F,T); phir(:,UN_nmf==1) = angle(X(:,UN_nmf==1));
        Xc_nmf = Vapprox .* exp(1i * phir);
        xc_nmf = real(iSTFT(Xc_nmf,Nfft,hop,Nw,wtype));
        
        % Phase unwrapping - oracle 
        phi = phase_unwrapping(X,hop,delta_oracle);
        XPU_oracle = V.* exp(1i*phi);
        xPU_oracle = real(iSTFT(XPU_oracle,Nfft,hop,Nw,wtype));
        
        % Phase unwrapping - nmf 
        phi = phase_unwrapping(Xapprox,hop,delta_nmf); 
        XPU_nmf = Vapprox.* exp(1i*phi);
        xPU_nmf = real(iSTFT(XPU_nmf,Nfft,hop,Nw,wtype));
        
        % Griffin Lim - oracle 
        [xGL_oracle,XGL_oracle] = griffin_lim(XC_oracle,Nw,hop,iter_GL,delta_oracle);
        
        % Griffin Lim - nmf
        [xGL_nmf,XGL_nmf] = griffin_lim(Xc_nmf,Nw,hop,iter_GL,delta_nmf);
        
        
        %%% Record audio files and compute the score
        
        % Score
        SDR(nw,1,it) = GetSDR(xGL_oracle ,x'); icons(nw,1,it) = norm(XGL_oracle-STFT(xGL_oracle,Nfft,hop,Nw,wtype))^2;
        SDR(nw,2,it) = GetSDR(xGL_nmf ,x'); icons(nw,2,it) = norm(XGL_nmf-STFT(xGL_nmf,Nfft,hop,Nw,wtype))^2;
        SDR(nw,3,it) = GetSDR(xPU_oracle',x'); icons(nw,3,it) = norm(XPU_oracle-STFT(xPU_oracle,Nfft,hop,Nw,wtype))^2;
        SDR(nw,4,it) = GetSDR(xPU_nmf',x'); icons(nw,4,it) = norm(XPU_nmf-STFT(xPU_nmf,Nfft,hop,Nw,wtype))^2;
        SDR(nw,5,it) = GetSDR(xc_oracle ,x'); icons(nw,5,it) = norm(XC_oracle-STFT(xc_oracle,Nfft,hop,Nw,wtype))^2;
        SDR(nw,6,it) = GetSDR(xc_nmf ,x'); icons(nw,6,it) = norm(Xc_nmf-STFT(xc_nmf,Nfft,hop,Nw,wtype))^2;
        
        % Record
        audiowrite(strcat(direc,'_',int2str(it),'_orig.ogg'),x,Fs);
        audiowrite(strcat(direc,'_',int2str(it),'_',int2str(Nw),'_oracle_corr.ogg'),xc_oracle,Fs);
        audiowrite(strcat(direc,'_',int2str(it),'_',int2str(Nw),'_oracle_GL.ogg'),xGL_oracle,Fs);
        audiowrite(strcat(direc,'_',int2str(it),'_',int2str(Nw),'_oracle_PU.ogg'),xPU_oracle,Fs);
        audiowrite(strcat(direc,'_',int2str(it),'_',int2str(Nw),'_nmf_corr.ogg'),xc_nmf,Fs);
        audiowrite(strcat(direc,'_',int2str(it),'_',int2str(Nw),'_nmf_GL.ogg'),xGL_nmf,Fs);
        audiowrite(strcat(direc,'_',int2str(it),'_',int2str(Nw),'_nmf_PU.ogg'),xPU_nmf,Fs);
    end
end

% Save scores
save(metrics_path,'gl_vs_pu.mat','SDR','icons');

% Plot results
sdr = mean(SDR,3); ic = mean(icons,3);
figure;
subplot(1,2,1); plot(1:LNw,sdr(:,1),'b-+',1:LNw,sdr(:,3),'r-*',1:LNw,sdr(:,5),'k-o',1:LNw,sdr(:,2),'b--+',1:LNw,sdr(:,4),'r--*',1:LNw,sdr(:,6),'k--o'); set(gca,'XtickLabel',wlength); xlabel('$N_w$ (samples)','fontsize',15,'interpreter','latex'); ylabel('SDR (dB)','fontsize',15);
subplot(1,2,2); semilogy(1:LNw,ic(:,1),'b-+',1:LNw,ic(:,3),'r-*',1:LNw,ic(:,5),'k-o',1:LNw,ic(:,2),'b--+',1:LNw,ic(:,4),'r--*',1:LNw,ic(:,6),'k--o'); set(gca,'XtickLabel',wlength); xlabel('$N_w$ (samples)','fontsize',15,'interpreter','latex');ylabel('Inconsistency','fontsize',15);
hl= legend('GL','PU','corrupted'); set(hl,'fontsize',13);