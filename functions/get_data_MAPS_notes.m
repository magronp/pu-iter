function [sm,x,Sm,X,ts,freq] = get_data_MAPS_notes(source_type,Fs,Nfft,Nw,hop,midi_pitch,wtype)

if nargin<7
    wtype = 'hann';
end

if nargin<6
    midi_pitch = randi(81,1,1)+20;
end

%%% Get data path
piano_path = 'datasets/MAPS/notes/piano';
switch source_type
    case 'rand'
        pitch2 = randi(88,1)+20;
        while pitch2 == midi_pitch
            pitch2 = randi(88,1)+20;
        end;
        source1_path = strcat(piano_path,int2str(midi_pitch),'_F.wav');
        source2_path = strcat(piano_path,int2str(pitch2),'_F.wav');
   case 'fifth'
       source1_path = strcat(piano_path,int2str(midi_pitch),'_F.wav');
       source2_path = strcat(piano_path,int2str(midi_pitch+7),'_F.wav');
end


%%% Read audio data

[s1,Fs_old] = audioread(source1_path); s1=s1';
s2 = audioread(source2_path); s2=s2';
s(1,:) = resample(s1,Fs,Fs_old); s(2,:) = resample(s2,Fs,Fs_old);
s = s(:,1:floor(size(s,2)/2));


%%% Build the mixture from the isolated sources
l = floor(size(s,2)); zz =zeros(1,l);
sm(1,:) = [s(1,:) zz  s(1,:) ];
sm(2,:) = [zz s(2,:)  s(2,:) ];

%%% Sources STFT
Sm = STFT(sm,Nfft,hop,Nw,wtype);
sm_aux = real(iSTFT(Sm,Nfft,hop,Nw,wtype));
sm = sm_aux;

%%% Mixture
X = sum(Sm,3); x = sum(sm,1);

[T] = size(X,2);
ts = (0:T-1)*hop / Fs;
freq = (1:Nfft/2)*Fs/Nfft;

end
