clear all; clc;

FileIndex = 3;

egm = load(['egm' num2str(FileIndex) '.csv']);

nde = size(egm,2); %256
N = 10173; %2035
m = 2; r = 0.1; % ApEn parameters
Fs = 2034.5; % sampling rate
[B A] = butter(3, [30/(Fs/2.0) 500/(Fs/2.0)]);


egm_filtered = zeros(N, nde);
for nodeId = 1:nde
    egm_filtered(:,nodeId) = filtfilt(B, A, squeeze(egm(:,nodeId)));
end


matlabpool open 4;
ApEn = zeros(1 ,nde);
parfor nodeId = 1:nde
    EGM = egm_filtered(:,nodeId);
    
    phi = zeros(1, 2);
    for j = 1:2
        mm = m+j-1;
        data = zeros(mm, N-mm+1);
        C = zeros(1, N-mm+1);
        
        for i = 1:mm
            data(i,:) = EGM(i:N-mm+i);
        end
        
        for i = 1:N-mm+1
            tmp = abs(data - repmat(data(:,i), 1, N-mm+1));
            C(i) = sum(~any((tmp > r), 1)) / (N-mm+1);
        end
        
        phi(j) = sum(log(C)) / (N-mm+1);
    end
    
    ApEn(nodeId) = phi(1) - phi(2);
end
matlabpool close;


fid = fopen(['ApEn' num2str(FileIndex) '.txt'], 'w');
fprintf(fid, '%f\n', ApEn);
fclose(fid);

f_=442; d_=1; fs_=44100; sound(sin(2*pi*f_*(1:d_*fs_)/fs_),fs_);