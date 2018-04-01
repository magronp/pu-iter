% Source separation procedure using the phase unwrapping technique
% initialization.
%
% Ref:
% "Model-based STFT phase recovery for audio source separation",
% Paul Magron, Roland Badeau and Bertrand David, TASLP, 2018
%
% Inputs:
%     X : F*T mixture
%     Xe : F*T*K sources (with amplitude and onset phases)
%     UN : K*T onset indicatrix function
%     hop : STFT overlap (in samples)
%     Nit : number of iterations
%     tplot : time frame where the mixing error is stored
%     unwr : =1 : initialisation with unwrapping ; =0 random
%
% Outputs:
%     Ye : estimated components
%     errormix : F*Nit matrix error in time frame 'tplot'


function [Ye,errormix] = phase_unwrap_ssep(X,Xe,UN,hop,Nit,tplot,unwr)

% Default values
if nargin<7, unwr=1; end
if nargin<6, tplot=0; end

% Parameters
[F,T,K] = size(Xe);
Nfft = (F-1)*2;
Ct = 2*pi*hop/Nfft;

% Weights: Wiener gain
lambda = abs(Xe).^2./(repmat(sum(abs(Xe).^2,3),[1 1 K])+eps);

% Initial values
Ye=Xe;
errormix = zeros(F,Nit+1);

% Loop over time frames
for t=2:T
    
        % Initialization
        for k=1:K
            if (UN(k,t)==0) % non-onset frame for source k
                if unwr
                    f_inf = freq_influence(abs(Ye(:,t,k))+eps)-1;
                    phiaux = angle(Ye(:,t-1,k))+Ct*f_inf;
                    Ye(:,t,k) = abs(Ye(:,t,k)) .* exp(1i * phiaux);
                else
                    Ye(:,t,k) = abs(Ye(:,t,k)) .* exp(1i * 2*pi*rand(F,1));
                end
            end
        end
        
        Yaux = squeeze(Ye(:,t,:));
        lambda_aux = squeeze(lambda(:,t,:));
            
        % Iiterative procedure
        E = X(:,t) - sum(Yaux,2);
                if t==tplot
                    errormix(:,1)=abs(E).^2;
                end
        v = abs(Yaux);
        for iter =1:Nit
            Yaux = Yaux + repmat(E,[1 K]).*lambda_aux;
            Yaux = Yaux ./ (abs(Yaux)+eps) .* v;
            E = X(:,t) - sum(Yaux,2);
                if t==tplot
                    errormix(:,iter+1)=abs(E).^2;
                end
        end
        
        % Update the sources in frame t
        Ye(:,t,:) = Yaux;
        
end

end


% PU: frequencies and regions of influence
function [f_inf,f_centr,f_harm] = freq_influence(v)

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
    f_inf = (1:length(v))'-1;
end

end

% Quadratic Interpolated FFT
function [p,b,a] = qint(ym1,y0,yp1)

p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1));
b = y0 - 0.25*(ym1-yp1)*p;
a = 0.5*(ym1 - 2*y0 + yp1);

end