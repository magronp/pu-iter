function [Xe] = consistent_wiener(X,V2,gamma,Nfft,Nw,hop,wtype)

% V2 : variances |STFT|^2 
% gamma : consistency weight
% Nfft : number of FFT points
% Nw : STFT window length
% hop : hop size (in samples)
% wtype = window type (Hann, Hamming...)

if nargin<7
    wtype = 'hann';
end


[F,T,K] = size(V2);
Xe = zeros(F,T,K);
mixup = X;

for k=1:K-1
    Vsources = V2(:,:,k:end);
    VS = Vsources(:,:,1);
    VN = sum(Vsources,3)-VS;
    
    aux = cons_wiener_penalty(mixup,VS+eps,VN+eps,gamma,Nfft,Nw,hop,wtype);
    Xe(:,:,k) = aux;
    mixup = mixup - aux;
end
Xe(:,:,K) = X-sum(Xe(:,:,1:K-1),3);

end



function [SE,iter] = cons_wiener_penalty(X,VS,VN,gamma,Nfft,Nw,hop,wtype)

% CONSWIENER_PENALTY Wiener filtering with STFT consistency penalty
%
% SE=conswiener_penalty(X,VS,VN,gamma,nsampl,hop)
%
% Inputs:
% X: nbin x nfram STFT-domain mixture signal
% VS: nbin x nfram STFT-domain source variance
% VN: nbin x nfram STFT-domain noise variance
% gamma: weight of the STFT consistency penalty
% hop: STFT hopsize
%
% Outputs:
% SE: nbin x nfram STFT-domain estimated source signal
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2013 Emmanuel Vincent and Jonathan Le Roux
%
% Patent Applied For in Japan: JP2011170190 
%
% Commercial use of this software may be subject to limitations.
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

[nbin,nfram]=size(X);

%%% Unconstrained Wiener filter %%%
SE=VS./(VS+VN).*X;
Lambda=1./(VS)+1./(VN);

%%% Conjugate gradient %%%
invM=1./(Lambda+gamma*(Nw-hop)/Nw);
wei=repmat([1; 2*ones(nbin-2,1); 1],[1 nfram]);
se=real(iSTFT(SE,Nfft,hop,Nw,wtype));
FSE=SE-STFT(se,Nfft,hop,Nw,wtype);
r=-gamma*FSE;
z=invM.*r;
P=z;
rsold=real(sum(sum(wei.*conj(r).*z)));
iter=0;
converged=false;
while ~converged,
    iter=iter+1;
    p=real(iSTFT(P,Nfft,hop,Nw,wtype));
    FP=P-STFT(p,Nfft,hop,Nw,wtype);
    AP=Lambda.*P+gamma*FP;
    alpha=rsold/real(sum(sum(wei.*conj(P).*AP))+realmin);
    SE=SE+alpha*P;
    converged=(sum(sum(alpha^2*real(P.*conj(P)))) < 1e-6*sum(sum(real(SE.*conj(SE)))));
    r=r-alpha*AP;
    z=invM.*r;
    rsnew=real(sum(sum(wei.*conj(r).*z)));
    beta=rsnew/(rsold+realmin);
    P=z+beta*P;
    rsold=rsnew;
end

end