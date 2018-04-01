clear all; close all; clc;
test_or_dev = 'Test';
set_settings_puiter;
selec = 'blind';

direc = strcat(audio_path,'separation_',selec,'/');

% PEASS options
options.segmentationFactor = 1;
options.destDir = direc;
score = zeros(3,4,J,Nsongs);

for it=1:Nsongs

    % Original Files path
    originalFiles = cell(J,1);
    for j=1:J
        originalFiles{j} = strcat(direc,'song',int2str(it),'_source',int2str(j),'_orig.wav');
    end

    % PEASS
    for j=1:J
        clc; fprintf('data %d / %d \n source %d / %d',it,Ndata,j,J);
        tic;
        % Wiener
        est_cnmf = strcat(direc,'song',int2str(it),'_source',int2str(j),'_',algos{1},'.wav');
        res = PEASS_ObjectiveMeasure(originalFiles,est_cnmf,options);
        score(1,:,j,it) = [res.OPS res.TPS res.IPS res.APS];

        % Consistent Wiener
        est_cnmf = strcat(direc,'song',int2str(it),'_source',int2str(j),'_',algos{2},'.wav');
        res = PEASS_ObjectiveMeasure(originalFiles,est_cnmf,options);
        score(2,:,j,it) = [res.OPS res.TPS res.IPS res.APS];
 
         % PU-iter
        est_cnmf = strcat(direc,'song',int2str(it),'_source',int2str(j),'_',algos{3},'.wav');
        res = PEASS_ObjectiveMeasure(originalFiles,est_cnmf,options);
        score(3,:,j,it) = [res.OPS res.TPS res.IPS res.APS]; 
        
        originalFiles = circshift(originalFiles,1);
    end
        
    delete(strcat(direc,'*eArtif.wav'));
    delete(strcat(direc,'*eInterf.wav'));
    delete(strcat(direc,'*eTarget.wav'));
    delete(strcat(direc,'*true.wav'));

end

save(strcat(metrics_path,'ssep_peass_',selec,'.mat'),'score');