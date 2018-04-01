% Set the settings used in the experiments for PU-Iter

Fs = 44100;

% Data
Nsongs = 50;
datavec = 1:Nsongs;
switch test_or_dev
    case 'Dev'
        dataNaN=[1 18 29 38 49];
    case 'Test'
        dataNaN=[6 34 36 40];
end
datavec(dataNaN)=[];
datavec = datavec(1:min(Nsongs,length(datavec)));
Nsongs = length(datavec);

J = 4;
L = 441344;  %song length after stft+istft (samples)

% STFT parameters
Nfft = 4096;            % number of FFT points
Nw = 4096;              % window length
hop = Nw/4;             % hop size
wtype = 'hann';         % window type

% Paths
dataset_path = 'datasets/DSD100/';
audio_path = 'audio_files/';
metrics_path = 'metrics/';

% Algorithms
algos = {'Wiener','ConsW','PU-Iter'};
Nalgo = length(algos);

% NMF
K = 50;                 % rank of the NMF per source
Ktot = K*J;             % total rank in the NMF model
iter_nmf = 200;         % number of iterations for the NMFs

% Consistent Wiener
gamma_wc = 4;

% PU-Iter
iter_puiter = 50;       % number of iterations for PU-Iter
