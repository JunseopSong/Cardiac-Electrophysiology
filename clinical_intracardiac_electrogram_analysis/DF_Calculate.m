clc; clear all;

%% Setting
Fs = 2034.5; % sampling rate
time = 5; % recording duration
numData = 256; % number of electrograms
BiEgmData = load('egm1.csv'); % input file
fid = fopen('DF1.txt', 'w'); % output file

%%
DF = zeros(numData, 1);
lenEgm = round(Fs*time);
nfft = 2^nextpow2(lenEgm);
[B A] = butter(3, [1/(Fs/2.0) 20/(Fs/2.0)]);
hannWindow = hanning(lenEgm);
timeVec = (0:lenEgm-1)/Fs*1000.0;
frequencyVec = Fs/2.0*linspace(0, 1, nfft/2+1);
freqIdxStart = find(frequencyVec>=3.0, 1, 'first');
freqIdxEnd = find(frequencyVec<=15.0, 1, 'last');

%%
for nodeId = 1:numData
    egm_raw = squeeze(BiEgmData(1:lenEgm, nodeId));
    
    egm_detrended = detrend(egm_raw, 'constant'); % baseline correction
    egm_rectified = abs(egm_detrended);
    egm_filtered = filtfilt(B, A, egm_rectified);
    
    % subplot(4,1,1); plot(timeVec, egm_raw); title('Bipolar electrogram raw data');
    % subplot(4,1,2); plot(timeVec, egm_detrended); title('Baseline drifting corrected');
    % subplot(4,1,3); plot(timeVec, egm_rectified); title('Rectified');
    % subplot(4,1,4); plot(timeVec, egm_filtered); title('1-20Hz third-order butterworth filtered');
    % xlabel('Time (ms)');
    
    signal_FT = fft(hannWindow.*egm_filtered, nfft)/lenEgm;
    spectrum_egm =  2.0*abs(signal_FT(1:nfft/2+1))/norm(hannWindow);
    
    % plot(frequencyVec, spectrum_egm); xlim([0 40]); xlabel('Frequency (Hz)');
    
    maxAmp = max(spectrum_egm(freqIdxStart:freqIdxEnd));
    DF(nodeId) = frequencyVec(find(spectrum_egm(freqIdxStart:freqIdxEnd) >= maxAmp) + freqIdxStart - 1);
    
    disp(nodeId);
end

fprintf(fid, '%f\n', DF);
fclose(fid);