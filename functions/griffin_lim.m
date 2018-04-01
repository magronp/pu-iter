%  Phase reconstruction algorithm based on [Griffin and Lim 1984]
% 
% Inputs :
%     X : complex matrix data to be modified
%     w : analysis and synthesis window
%     hop : temporal shift between 2 frames (number of samples)
%     n_iter : number of iterations
%     delta : F*T binary mask where GL updates are performed
%     wtype : (string) STFT window type
% 
% Output :
%     y : estimated time signal
%     Y : estimated STFT
%     err : vector of the square error between V and |X| over iteration
%     icons : inconsistency over iterations

function [y,Y,err,icons] = griffin_lim(X,Nw,hop,n_iter,delta,wtype)

if nargin<6
    wtype = 'hann';
end

if nargin<5
    delta = ones(size(X));
end

Nfft = (size(X,1)-1)*2;
V = abs(X);

% Erros
err = zeros(1,n_iter);
icons = zeros(1,n_iter);

% Initialization
Y = X; Yold=Y;
y = real(iSTFT(Y, Nfft, hop, Nw, wtype));

% G&L algorithm
for it=1:n_iter
    
    % STFT from y
    Ynew = STFT(y, Nfft,  hop, Nw, wtype);
    
    modY = abs(Ynew)+eps;
    err(it) = norm(V-modY);
    icons(it) = norm(Yold-Ynew)^2;
    
    % Target magnitude value
    Y = V.*Ynew./modY;
    Y(delta==0) = X(delta==0);
    
    y = real(iSTFT(Y, Nfft, hop, Nw, wtype));
    Yold = Ynew;

end

end