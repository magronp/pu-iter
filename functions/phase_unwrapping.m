% Linear phase unwrapping algorithm
% Paul Magron, Feb 2015
%
% Inputs:
%     X : F*T complex matrix (random or initial phase if any, whose position is given by the mask delta)
%     hop : STFT overlap (in samples)
%     delta : index matrices (the bin (f,t) is reconstructed if delta(f,t)=1)
%
% Outputs:
%     phi : reconstructed phase
%     f_inf : frequencies

function [phi,f_inf] = phase_unwrapping(X,hop,delta)

V = abs(X);
[F,T] = size(V);

% if no mask provided, reconstruct phase everywhere
if nargin<3
    delta = ones(F,T);
end

% Get frequencies and decomposition in regions of influence
f_inf = get_frequencies_qifft(V);
Nfft = 2*(F-1);
Ct = 2*pi*hop/Nfft;

% Phase unwrapping
phi = angle(X);
for t=2:T
    phi_pu = phi(:,t-1)+Ct*f_inf(:,t);
    % Only use PU if delta=1 (else, we keep the phase in X)
    phi(delta(:,t)==1,t) = phi_pu(delta(:,t)==1);
end
phi = angle(exp(1i*phi));

end

% Frequencies and decomposition in regions of influence
function [f_inf,f_centr,f_harm] = get_frequencies_qifft(V)

[F,T] = size(V);

% Instantaneous frequencies
f_inf = zeros(F,T);
f_centr = zeros(F,T);
f_harm = cell(1,T);

for t=1:T
    
    v = V(:,t)';

    % Magnitude peaks
    [~,centr] = findpeaks(max(v,10^-6),'MINPEAKHEIGHT',max(v)*10^-2);

    Nfreq = length(centr);
    harm = zeros(1,Nfreq);

    if (Nfreq >0)
        % Quadratic interpolation of frequencies
        for ind = 1:Nfreq
            f = centr(ind);
            harm(ind) = qint((log(v(f-1))),(log(v(f))),(log(v(f+1))))+f;
        end

        % Frequencies in Regions of influence
        inf = zeros(length(v),1);
        deb = 1;
        index_lim = zeros(1,Nfreq-1);

        for ind = 1:(Nfreq-1)
            f = centr(ind);
            fp = centr(ind+1);
            fin = floor((v(fp)*f+v(f)*fp)/(v(fp)+v(f)));
            inf(deb:fin) = harm(ind);
            deb = fin+1;
            index_lim(ind) = fin;
        end

        inf(deb:end) = harm(end);
        f_centr(centr(:),t) = 1;
        
    else
        inf = (1:length(v))'-1;
        f_centr(:,t) = 0;
    end

    f_inf(:,t) = inf-1;
    f_harm{t} = harm;
end
end

% QIFFT
function [p,b,a] = qint(ym1,y0,yp1)

p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1));
b = y0 - 0.25*(ym1-yp1)*p;
a = 0.5*(ym1 - 2*y0 + yp1);

end
