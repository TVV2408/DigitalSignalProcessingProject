clear; clc; close all;

% ===== Paths =====
rootDir     = 'data';
cleanPath   = fullfile(rootDir, 'clean');
mixedPath   = fullfile(rootDir, 'mixed');
outPath     = fullfile('outputs', 'lms_filtered_fix');
if ~exist(outPath, 'dir'), mkdir(outPath); end

% ===== Options =====
Fs_target   = 16000;       % tần số lấy mẫu
M           = 128;         % số tap FIR / độ dài bộ lọc thích nghi
mu          = 0.002;       % tốc độ cập nhật trọng số / step size
useNLMS     = true;        % bật/tắt NLMS
alpha       = 0.995;       % hệ số rò rỉ (leaky)
wClip       = 10;          % chặn biên trọng số
antiClipPK  = 0.999;       
printLog    = true;

cleanList = dir(fullfile(cleanPath, '*.flac'));
mixedList = dir(fullfile(mixedPath, '*_snr*dB.wav'));
assert(~isempty(mixedList), 'Không thấy file .wav trong data/mixed/');
assert(~isempty(cleanList), 'Không thấy file .flac trong data/clean/');

cleanMap = containers.Map();
for k = 1:numel(cleanList)
    [~,nm,~] = fileparts(cleanList(k).name);
    cleanMap(lower(nm)) = fullfile(cleanList(k).folder, cleanList(k).name);
end

if printLog
    fprintf('BẮT ĐẦU LMS cho %d file mixed...\n', numel(mixedList));
    fprintf('%-32s %-8s %-8s %-8s %-8s\n', 'File', 'SNR_b', 'SNR_a', 'ΔSNR', 'RMS↓%');
end

for i = 1:numel(mixedList)
    mixFile = fullfile(mixedList(i).folder, mixedList(i).name);
    [x_noisy, Fs_m] = audioread(mixFile); x_noisy = mean(x_noisy,2);

    [~, base, ~] = fileparts(mixedList(i).name);    
    parts = split(string(base), '_');
    if numel(parts) < 3, continue; end
    cleanBase = char(lower(parts(1)));
    if ~isKey(cleanMap, cleanBase), continue; end
    cleanFile = cleanMap(cleanBase);

    [x_clean, Fs_c] = audioread(cleanFile); x_clean = mean(x_clean,2);

    if Fs_m ~= Fs_target, x_noisy = resample(x_noisy, Fs_target, Fs_m); end
    if Fs_c ~= Fs_target, x_clean = resample(x_clean, Fs_target, Fs_c); end
    Fs = Fs_target;

    L = min(length(x_noisy), length(x_clean));
    x_noisy = x_noisy(1:L);
    x_clean = x_clean(1:L);

    % ===== LMS / NLMS =====
    w = zeros(M,1);
    y = zeros(L,1);
    e = zeros(L,1);
    for n = M:L
        xvec = x_noisy(n:-1:n-M+1);
        y(n) = w.' * xvec;
        e(n) = x_clean(n) - y(n);
        if useNLMS
            denom = (xvec.'*xvec) + 1e-8;
            w = alpha*w + (mu/denom) * xvec * e(n);
        else
            w = alpha*w + mu * xvec * e(n);
        end
        if any(abs(w) > wClip)
            w = wClip * tanh(w / wClip);
        end
    end
    x_hat = e;        

    pk = max(abs(x_hat));
    if pk > antiClipPK, x_hat = antiClipPK * x_hat / pk; end

    e_in  = x_noisy - x_clean;
    e_out = x_hat   - x_clean;
    SNRb  = 10*log10(sum(x_clean.^2)/(sum(e_in.^2)+eps));
    SNRa  = 10*log10(sum(x_clean.^2)/(sum(e_out.^2)+eps));
    RMSb  = sqrt(mean(e_in.^2));
    RMSa  = sqrt(mean(e_out.^2));
    RMSred = (1 - RMSa/max(RMSb,eps))*100;

    % Ghi file
    outFile = fullfile(outPath, [base '_lms.wav']);
    audiowrite(outFile, x_hat, Fs);

    if printLog
        fprintf('%-32s %-8.2f %-8.2f %-8.2f %-8.2f\n', [base '.wav'], SNRb, SNRa, (SNRa-SNRb), RMSred);
    end
end

fprintf('--- Hoàn tất LMS filtering ---\n');
