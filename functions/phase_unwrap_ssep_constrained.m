% Source separation procedure using the phase unwrapping technique
% initialization and a relaxed phase constraint
% Paul Magron, May 2016
%
% Inputs:
%     X : F*T mixture
%     Xe : F*T*K sources (with amplitude and onset phases)
%     UN : K*T onset indicatrix function
%     hop : STFT overlap (in samples)
%     Nit : number of iterations
%     sigma = constraint weight parameter
%     unwr : =1 : initialisation with unwrapping ; =0 random
%     Yrand : initial random values
%
% Outputs:
%     Ye : estimated components


function [Ye,errormix] = phase_unwrap_ssep_constrained(X,Xe,UN,hop,Nit,sigma,unwr,Yrand)

% Parameters
[F,T,K] = size(Xe);
Nfft = (F-1)*2;
Ct = 2*pi*hop/Nfft;

% Weights
G = abs(Xe).^2./(repmat(sum(abs(Xe).^2,3),[1 1 K])+eps);

% Initial values
Ye=Xe;
Z= Xe;
errormix = zeros(F,Nit+1);

% Input default values
if nargin<8
    Yrand = abs(Ye) .* exp(1i * 2*pi*rand(size(Ye)));
end

if nargin<7
    unwr=1;
end


% Loop over time frames
for t=2:T   
        % Initialisation : if onset frame, do nothing, if not, unwrapping
        % or random values
        for k=1:K
            if (UN(k,t)==0) % onset frame for source k
                 f_inf = get_frequencies_qifft_frame(abs(Ye(:,t,k)));
                 phiaux = angle(Ye(:,t-1,k))+Ct*f_inf;
                 Z(:,t,k) = abs(Ye(:,t,k)) .* exp(1i * phiaux);
                if unwr
                    Ye(:,t,k) =Z(:,t,k);
                else
                    Ye(:,t,k) = Yrand(:,t,k);
                end
            end
        end
        
        % Estimation from the mixture
        Yaux = squeeze(Ye(:,t,:));
        Zaux = squeeze(Z(:,t,:));
        Gaux = squeeze(G(:,t,:));

        E = X(:,t) - sum(Yaux,2);
        for iter =1:Nit
            Yaux = Yaux + repmat(E,[1 K]).*Gaux;
            Yaux = (Yaux +  sigma .* Gaux .* Zaux  ) ./ (abs(Yaux  +  sigma.* Gaux .* Zaux  )+eps) .* abs(Yaux);
            E = X(:,t) - sum(Yaux,2);
        end
        Ye(:,t,:) = Yaux;
        
end

end

% PU: frequencies and regions of influence
function [f_inf,f_centr,f_harm] = get_frequencies_qifft_frame(v)

v = v(:)';

%Central peaks
%[~,f_centr] = findpeaks(v,'MINPEAKHEIGHT',0.01*max(v));
[~,f_centr] = findpeaks(max(v,10^-6),'MINPEAKHEIGHT',max(v)*10^-2);

Nfreq = length(f_centr);
f_harm = zeros(1,Nfreq);

if (Nfreq >0)
    % QIFFT
    for ind = 1:Nfreq
        f = f_centr(ind);
        f_harm(ind) = qint(log10(v(f-1)),log10(v(f)),log10(v(f+1)))+f;
    end

    % Regions of influence
    f_inf = zeros(length(v),1);
    deb = 1;
    index_lim = zeros(1,Nfreq-1);
    
    for ind = 1:(Nfreq-1)
        f = f_centr(ind);
        fp = f_centr(ind+1);
        fin = floor((v(fp)*f+v(f)*fp)/(v(fp)+v(f)));
        f_inf(deb:fin) = f_harm(ind);
        deb = fin+1;
        index_lim(ind) = fin;
    end

    f_inf(deb:end) = f_harm(end);
    
else
    f_inf = (1:length(v))';
end

%  Frequencies start from 0
f_inf = f_inf-1;

end

% Quadratic Interpolated FFT
function [p,b,a] = qint(ym1,y0,yp1)

p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1));
b = y0 - 0.25*(ym1-yp1)*p;
a = 0.5*(ym1 - 2*y0 + yp1);

end
