% Source separation benchmark on DSD100 database using several phase
% recovery methods

clc; clear all; close all;
test_or_dev = 'Test';
set_settings_puiter;
selec = 'blind';

SDR = zeros(J,Nalgo,Nsongs); SIR = zeros(J,Nalgo,Nsongs); SAR = zeros(J,Nalgo,Nsongs);
time_comput = zeros(Nsongs,Nalgo);

for it =1:Nsongs
   
    % Source generation
    clc; fprintf('Data %d / %d \n',it,Nsongs);
    num_piece = datavec(it);
    [sm,x,Sm,X,ts,freq] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
    [F,T,J] = size(Sm);
    V = abs(Sm);
    Xe = zeros(F,T,J,Nalgo);
    xe = zeros(J,length(x),Nalgo);
    
    % Magnitude estimation
    fprintf('Magnitude estimation... \n');
    switch selec
        case 'oracle'
            Vapprox = V;
        case 'informed'
            Wini=rand(F,K); Hini=rand(K,T);
            Vapprox = zeros(F,T,J);
            for inst=1:J
                [waux,haux] = NMF(abs(Sm(:,:,inst)),Wini,Hini,iter_nmf,1,0);
                Vapprox(:,:,inst) = waux*haux;
            end
            
        case 'blind'
            Wini=rand(F,Ktot); Hini=rand(Ktot,T);
            [W,H,cost] = NMF(abs(X),Wini,Hini,iter_nmf,1,0);
            genom=Clustering(W,H,J);
            Vapprox = zeros(F,T,J);
            for j=1:J
                Vapprox(:,:,j) = W(:,genom==j)*H(genom==j,:);
            end
    end
    
    V2 = Vapprox.^2;
    Sm_approx = Vapprox .* exp(1i*repmat(angle(X),[1 1 J]));
    
    % Wiener filtering
    tic;
    Xe(:,:,:,1) = V2 .* repmat(X ./ (sum(V2,3)+eps),[1 1 J]);
    time_comput(it,1) = toc;
    
    % Consistent Wiener filtering
    fprintf('Consistent Wiener filtering... \n');
    tic;
    Xe(:,:,:,2) = consistent_wiener(X,V2,gamma_wc,Nfft,Nw,hop,wtype);
    time_comput(it,2) = toc;
    
    % PU-Iter
    fprintf('Iter... \n');
    UN = detect_onset_frames(Vapprox,Fs,hann(Nw),hop);
    tic;
    Xe(:,:,:,3) = phase_unwrap_ssep(X,Sm_approx,UN,hop,iter_puiter);
    time_comput(it,3) = toc;
    
    %%% Synthesis and record
    
    % iSTFT
    for al=1:Nalgo
        xe(:,:,al) = real(iSTFT(squeeze(Xe(:,:,:,al)),Nfft,hop,Nw,wtype));
    end
        
    % Score
    for al=1:Nalgo
        [sdr,sir,sar,perm] = GetSDR(squeeze(xe(:,:,al)),sm);
        SDR(:,al,it) = sdr; SIR(:,al,it) = sir; SAR(:,al,it) = sar;
    end
    xe_perm = xe(perm,:,:);
   
    % Record
    fprintf('Record... \n');
    audiowrite(strcat(audio_path,'separation_',selec,'/song',int2str(it),'.wav'),x,Fs);
    for k=1:J
        audiowrite(strcat(audio_path,'separation_',selec,'/song',int2str(it),'_source',int2str(k),'_orig.wav'),sm(k,:),Fs);
        for al=1:Nalgo
            audiowrite(strcat(audio_path,'separation_',selec,'/song',int2str(it),'_source',int2str(k),'_',algos{al},'.wav'),squeeze(xe_perm(k,:,al)),Fs);
        end
    end

    
end

% Save the BSS eval scores and computation time
save(strcat('phase unwrapping source sep/benchmark/score_bss_',selec,'.mat'),'SDR','SIR','SAR','time_comput');
