% Onset frames detection from a spectrogram
%
% Inputs :
%     V : spectrogram (F*T matrix)
%     Fs : sample rate
%     w : analysis window
%     hop : hop size
%
% Output :
%     UN : K*T matrix =1 if onset and 0 else

function [UN] = detect_onset_frames(V,Fs,w,hop)

[F,T,K] = size(V);
UN = zeros(K,T);

% If source type unspecified, use Tempogram toolbox
Nfft = 2*(F-1);
parameter.win_len = Nfft;
parameter.fs = Fs;
parameter.stepsize = hop;
parameter.StftWindow = w;

for k=1:K
    [v0,~] = audio_to_noveltyCurve(V(:,:,k), parameter);
    [~,ind] = findpeaks(v0,'MINPEAKHEIGHT',max(v0)*0.1);
    t0 = ind'; t0 = t0(1:end);
    UN(k,t0)=1;
end

end


function [noveltyCurve,featureRate] = audio_to_noveltyCurve(V, parameter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a MODIFIED version of the original method by Grosche and Muller :
% http://www.mpi-inf.mpg.de/resources/MIR/tempogramtoolbox/
%
% Description:
% Computes a novelty curve (onset detection function) 
% for the input audio spectrogram. This implementation is a variant
% of the widely used spectral flux method with additional 
% bandwise processing and a logarithmic intensity compression.
% This particularly addresses music with weak onset information 
% (e.g., exhibiting string instruments.)
%
% Input:
%       V : spectrogram
%       fs : sampling rate of the audio (Hz)
%       parameter (optional): parameter struct with fields
%               .logCompression : enable/disable log compression
%               .compressionC : constant for log compression
%               .win_len : window length for STFT (in samples)
%               .stepsize : stepsize for the STFT
%               .resampleFeatureRate : feature rate of the resulting
%               novelty curve (resampled, independent of stepsize)
%
% Output:
%       noveltyCurve : the novelty curve
%       featureRate : feature rate of the novelty curve (Hz)
%
% License:
%     This file is part of 'Tempogram Toolbox'.
% 
%     'Tempogram Toolbox' is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 2 of the License, or
%     (at your option) any later version.
% 
%     'Tempogram Toolbox' is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with 'Tempogram Toolbox'. If not, see
%     <http://www.gnu.org/licenses/>.
%
% 
% Reference:
%   Peter Grosche and Meinard M�ller
%   Extracting Predominant Local Pulse Information from Music Recordings 
%   IEEE Transactions on Audio, Speech, and Language Processing, 19(6), 1688-1701, 2011.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    parameter = [];
end

if isfield(parameter,'fs')==0
    parameter.fs = 22050;
end
if isfield(parameter,'win_len')==0
    parameter.win_len = 1024*parameter.fs/22050;
end
if isfield(parameter,'stepsize')==0
    parameter.stepsize = 512*parameter.fs/22050;
end
if isfield(parameter,'compressionC')==0
    parameter.compressionC = 1000;
end
if isfield(parameter,'logCompression')==0
    parameter.logCompression = true;
end
if isfield(parameter,'resampleFeatureRate')==0
    parameter.resampleFeatureRate = 200;
end

fs = parameter.fs;

myhann = @(n)  0.5-0.5*cos(2*pi*((0:n-1)'/(n-1)));
%myhann = @(Nw)  perfect_reco_window(Nw,Nw/4);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameter.returnMagSpec = 1;
%[specData,featureRate] = audio_to_spectrogram_via_STFT(f_audio,parameter);
hop = parameter.stepsize;
featureRate = fs/hop;
parameter.featureRate = featureRate;
specData = V;

% normalize and convert to dB
specData = specData./max(max(specData));
thresh = -74; % dB
thresh = 10^(thresh./20);
specData = (max(specData,thresh));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bandwise processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bands = [0 500; 500 1250; 1250 3125; 3125 7812.5; 7812.5 floor(parameter.fs./2)]; %hz
compressionC = parameter.compressionC;

bandNoveltyCurves = zeros(size(bands,1),size(specData,2));

for band = 1:size(bands,1)

    bins = round(bands(band,:)./ (parameter.fs./parameter.win_len));
    bins = max(1,bins);
    bins = min(round(parameter.win_len./2)+1,bins);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % band novelty curve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    bandData = specData(bins(1):bins(2),:);
    if parameter.logCompression && parameter.compressionC>0
        bandData = log(1 + bandData.*compressionC)/(log(1+compressionC)); 
    end
    
    % smoothed differentiator
    diff_len = 0.3;%sec
    diff_len = max(ceil(diff_len*parameter.fs/parameter.stepsize),5);
    diff_len = 2*round(diff_len./2)+1;
    diff_filter = myhann(diff_len).*[-1*ones(floor(diff_len/2),1); 0;ones(floor(diff_len/2),1) ];
    diff_filter = diff_filter(:)';    
    bandDiff = filter2(diff_filter, [repmat(bandData(:,1),1,floor(diff_len/2)),bandData,repmat(bandData(:,end),1,floor(diff_len/2))]);
    bandDiff = bandDiff.*(bandDiff>0);
    bandDiff = bandDiff(:,floor(diff_len/2):end-floor(diff_len/2)-1);

    
    % normalize band
    norm_len = 5;%sec
    norm_len = max(ceil(norm_len*parameter.fs/parameter.stepsize),3);
    norm_filter = myhann(norm_len);
    norm_filter = norm_filter(:)';
    norm_curve = filter2(norm_filter./sum(norm_filter),sum(bandData));
    % boundary correction
%     norm_filter_sum = (sum(norm_filter)-cumsum(norm_filter))./sum(norm_filter);
%     norm_curve(1:floor(norm_len/2)) = norm_curve(1:floor(norm_len/2))./ fliplr(norm_filter_sum(1:floor(norm_len/2)));
%     norm_curve(end-floor(norm_len/2)+1:end) = norm_curve(end-floor(norm_len/2)+1:end)./ norm_filter_sum(1:floor(norm_len/2));
    bandDiff = bsxfun(@rdivide,bandDiff,norm_curve);

    % compute novelty curve of band 
    noveltyCurve = sum( bandDiff);
    bandNoveltyCurves(band,:) = noveltyCurve;

end

% figure; 
% for band = 1:5
%     subplot(5,1,6-band);plot(bandNoveltyCurves(band,:));
%     ylim([0 1])
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% summary novelty curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noveltyCurve = mean(bandNoveltyCurves);


% resample curve
% if (parameter.resampleFeatureRate >0 && parameter.resampleFeatureRate ~= parameter.featureRate)
%     [noveltyCurve,featureRate] = resample_noveltyCurve(noveltyCurve,parameter);
%     parameter.featureRate = featureRate;
% 
% end


% local average subtraction
[noveltyCurve] = novelty_smoothedSubtraction(noveltyCurve,parameter);



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [noveltySub,local_average] = novelty_smoothedSubtraction(noveltyCurve,parameter)

    myhann = @(n)  0.5-0.5*cos(2*pi*((0:n-1)'/(n-1)));
    smooth_len = 1.5;%sec
    smooth_len = max(ceil(smooth_len*parameter.fs/parameter.stepsize),3);
    smooth_filter = myhann(smooth_len);
    smooth_filter = smooth_filter(:)';
    local_average = filter2(smooth_filter./sum(smooth_filter),noveltyCurve);
    
    noveltySub = (noveltyCurve-local_average);
    noveltySub = (noveltySub>0).*noveltySub;
end

function [noveltyCurve,featureRate] = resample_noveltyCurve(noveltyCurve,parameter)

    p = round(1000*parameter.resampleFeatureRate./parameter.featureRate);
    noveltyCurve = resample(noveltyCurve,p,1000,10);
    featureRate = parameter.featureRate*p/1000;
end