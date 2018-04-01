clc; clear all; close all;
test_or_dev = 'Test';
set_settings_puiter;

% Load data
[sm,x,Sm,X,ts,freq] = get_data_MAPS_notes('fifth',Fs,Nfft,Nw,hop,60);
[F,T,K] = size(Sm);
V = abs(Sm);

Sm_approx = V .* exp(1i*repmat(angle(X),[1 1 K]));

% Onset detection
UN = detect_onset_frames(abs(Sm_approx),Fs,hann(Nw),hop);

% Iter - random and PU init
Nit = 50;
Xr = phase_unwrap_ssep(X,Sm_approx,UN,hop,Nit,100,0);
Xu = phase_unwrap_ssep(X,Sm_approx,UN,hop,Nit,100,1);


% Wiener and Consistent Wiener
V2 = V.^2;
Xw = V2 .* repmat(X ./ (sum(V2,3)+eps),[1 1 K]);
Xc = consistent_wiener(X,V2,4,Nfft,Nw,hop);

% Synthesis
xm=real(iSTFT(Sm_approx,Nfft,hop,Nw,wtype));
xr=real(iSTFT(Xr,Nfft,hop,Nw,wtype));
xu=real(iSTFT(Xu,Nfft,hop,Nw,wtype));
xw=real(iSTFT(Xw,Nfft,hop,Nw,wtype));
xc=real(iSTFT(Xc,Nfft,hop,Nw,wtype));


% STFTs with hop = 1
Xx = STFT(sm(1,:),Nfft,1,Nw,wtype);
Xmm = STFT(xm(1,:),Nfft,1,Nw,wtype);
Xrr = STFT(xr(1,:),Nfft,1,Nw,wtype);
Xuu = STFT(xu(1,:),Nfft,1,Nw,wtype);
Xww = STFT(xw(1,:),Nfft,1,Nw,wtype);
Xcc = STFT(xc(1,:),Nfft,1,Nw,wtype);

% Real parts
ttdeb =110000; ttend =110250;
tss = (0:size(Xx,2)-1)*1 / Fs;
ttts = tss(ttdeb:ttend);
fplot = 74;
vo = real(Xx(fplot,ttdeb:ttend));
vm = real(squeeze(Xmm(fplot,ttdeb:ttend)));
vr = real(squeeze(Xrr(fplot,ttdeb:ttend)));
vu = real(squeeze(Xuu(fplot,ttdeb:ttend)));
vw = real(squeeze(Xww(fplot,ttdeb:ttend)));
vc = real(squeeze(Xcc(fplot,ttdeb:ttend)));

% Compare Iter initializations (rand/PU)
figure;
plot(ttts,vo,'k-o',ttts,vm,'bx-',ttts,vr,'m-+',ttts,vu,'r*-');
xlabel('Time (s)','FontSize',16);
ylabel('Real part','FontSize',16);
Vax=axis; axis([ttts(1) ttts(end) Vax(3) Vax(4)]);
hl=legend('Original','Mixture','Random','PU'); set(hl,'FontSize',14);


% Compare Iter and Wiener/Cons Wiener
figure;
plot(ttts,vo,'k+-',ttts,vw,'bo-',ttts,vc,'ms-',ttts,vu,'r*-');
xlabel('Time (s)','FontSize',16);
ylabel('Real part','FontSize',16);
Vax=axis; axis([ttts(1) ttts(end) Vax(3) Vax(4)]);
hl=legend('Original','Wiener','Cons-W','PU-Iter'); set(hl,'FontSize',14);
