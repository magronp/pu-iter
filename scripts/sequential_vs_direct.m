% Compare iterative algorithms using a direct PU on the whole matrix or a
% sequential approach
clc; clear all; close all;
test_or_dev = 'Test';
set_settings_puiter;

diff_algo = 2;
SDR = zeros(Nsongs,diff_algo); SIR = zeros(Nsongs,diff_algo); SAR = zeros(Nsongs,diff_algo);

for it =1:Nsongs
    
    % Source generation
    clc; fprintf('Data %d / %d \n',it,Nsongs);
    num_piece = it;
    [sm,x,Sm,X,ts,freq] = get_data_DSD(dataset_path,test_or_dev,num_piece,Fs,Nfft,Nw,hop);
    [F,T,J] = size(Sm);
    xe = zeros(J,length(x),diff_algo);
    
    % Onset detection
    UN = detect_onset_frames(abs(Sm),Fs,hann(Nw),hop);
    
    % Sequential algorithm
    tic;
    X_loop = phase_unwrap_ssep(X,Sm,UN,hop,iter_puiter);
    t_loop=toc;
    
    % Direct phase unwrapping + iter procedure (matricial)
    tic;
    X_direct = phase_unwrap_ssep(X,Sm,UN,hop,0,0);
    Xrepet = repmat(X,[1 1 J]);
    G = abs(Sm).^2 ./ (repmat( sum(abs(Sm).^2,3),[1 1 J])+eps);
    for iter=1:iter_puiter
        E = X - sum(X_direct,3);
        Y = X_direct + G .* repmat(E,[1 1 J]);
        X_direct =  Y./(abs(Y)+eps).*abs(Sm);
    end
    t_direct=toc;
    
    
    % Synthesis
    for j=1:J
        xe(j,:,1)=real(iSTFT(squeeze(X_loop(:,:,j)),Nfft,hop,Nw,wtype));
        xe(j,:,2)=real(iSTFT(squeeze(X_direct(:,:,j)),Nfft,hop,Nw,wtype));
    end
    
    % Score
    for al=1:diff_algo
        [sdr,sir,sar] = GetSDR(squeeze(xe(:,:,al)),sm);
        SDR(it,al) = mean(sdr); SIR(it,al) = mean(sir); SAR(it,al) = mean(sar);
    end
   
end

% Save score
save(strcat(metrics_path,'sequential_vs_direct.mat'),'SDR','SIR','SAR');
