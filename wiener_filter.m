clear; clc; close all;

rootDir    = 'data';
cleanPath  = fullfile(rootDir, 'clean');
mixedPath  = fullfile(rootDir, 'mixed');
outPath    = fullfile('outputs', 'wiener_filtered_fix');

if ~exist(outPath, 'dir'); mkdir(outPath); end

% Options
Fs_target  = 16000;               
win        = hamming(512, 'periodic'); 
hop        = 256;                 
nfft       = 512;

cleanList = dir(fullfile(cleanPath, '*.flac'));
mixedList = dir(fullfile(mixedPath, '*.wav'));
assert(~isempty(mixedList), 'Không thấy file .wav trong data/mixed/');
assert(~isempty(cleanList), 'Không thấy file .flac trong data/clean/');

cleanMap = containers.Map();
for k = 1:numel(cleanList)
    cleanMap(lower(erase(cleanList(k).name, '.flac'))) = fullfile(cleanList(k).folder, cleanList(k).name);
end

fprintf('Bắt đầu Wiener filtering cho %d file mixed...\n', numel(mixedList));

% ============
for i = 1:numel(mixedList)
    mixFile = fullfile(mixedList(i).folder, mixedList(i).name);
    [x_noisy, Fs_mix] = audioread(mixFile);
    x_noisy = mean(x_noisy, 2); 

    [~, baseName, ~] = fileparts(mixedList(i).name);
    parts = split(baseName, '_');
    if numel(parts) < 1
        fprintf('Bỏ qua file tên lạ: %s\n', mixedList(i).name);
        continue;
    end
    cleanBase = lower(string(parts{1})); 

    if ~isKey(cleanMap, cleanBase)
        fprintf('Không tìm thấy clean cho %s -> %s.flac\n', mixedList(i).name, cleanBase);
        continue;
    end
    cleanFile = cleanMap(cleanBase);
    [x_clean, Fs_clean] = audioread(cleanFile);
    x_clean = mean(x_clean, 2);

    if Fs_mix ~= Fs_target
        x_noisy = resample(x_noisy, Fs_target, Fs_mix);
    end
    if Fs_clean ~= Fs_target
        x_clean = resample(x_clean, Fs_target, Fs_clean);
    end
    Fs = Fs_target;

    L = min(length(x_noisy), length(x_clean));
    x_noisy = x_noisy(1:L);
    x_clean = x_clean(1:L);

    S_noisy = stft(x_noisy, Fs, 'Window', win, 'OverlapLength', hop, 'FFTLength', nfft);
    S_clean = stft(x_clean, Fs, 'Window', win, 'OverlapLength', hop, 'FFTLength', nfft);

    Pxx = abs(S_clean).^2;
    Pvv = abs(S_noisy - S_clean).^2;

    Hw  = Pxx ./ (Pxx + Pvv + eps);
    Y   = Hw .* S_noisy;

    x_hat = istft(Y, Fs, 'Window', win, 'OverlapLength', hop, 'FFTLength', nfft);
    
    Lhat = length(x_hat);
    if Lhat < L
        x_hat = [x_hat; zeros(L - Lhat, 1)];
    else
        x_hat = x_hat(1:L);
    end

    pk = max(abs(x_hat));
    if pk > 0.999
        x_hat = 0.999 * x_hat / pk;
    end

    e_in  = x_noisy - x_clean;                  
    e_out = x_hat   - x_clean;                   
    SNRin  = 10*log10(sum(x_clean.^2) / (sum(e_in.^2)  + eps));
    SNRout = 10*log10(sum(x_clean.^2) / (sum(e_out.^2) + eps));

    outFile = fullfile(outPath, sprintf('%s_wiener.wav', baseName));
    audiowrite(outFile, x_hat, Fs);

    fprintf('Done: %-40s | SNRin = %+5.1f dB -> SNRout = %+5.1f dB\n', ...
        [baseName '.wav'], SNRin, SNRout);
end

fprintf('--- Hoàn tất Wiener filtering ---\n');
