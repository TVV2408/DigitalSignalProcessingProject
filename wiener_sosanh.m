clear; clc; close all;

% Options
clean_path    = 'data/clean/';
mixed_path    = 'data/mixed/';
filtered_path = 'outputs/wiener_filtered_fix/';
Fs_target = 16000;      
doListen  = false;      
doFigure  = true;       
saveCSV   = true;       

mixed_files = dir(fullfile(mixed_path, '*_snr*dB.wav'));
clean_files = dir(fullfile(clean_path, '*.flac'));
if isempty(mixed_files), error('Không thấy file mixed trong %s', mixed_path); end
if isempty(clean_files), error('Không thấy file clean trong %s', clean_path); end

cleanMap = containers.Map();
for k = 1:numel(clean_files)
    [~,nm,~] = fileparts(clean_files(k).name);
    cleanMap(lower(nm)) = fullfile(clean_files(k).folder, clean_files(k).name);
end

fprintf('\n KẾT QUẢ ĐÁNH GIÁ WIENER FILTER (toàn bộ mixed)\n');
fprintf('----------------------------------------------------------------------------------------------\n');
fprintf('%-28s %-10s %-8s %-10s %-10s %-10s %-9s\n','File','Noise','SNRtag','SNR_b','SNR_a','DeltaSNR','RMS_down(%%)');
fprintf('----------------------------------------------------------------------------------------------\n');

M = {};  % {file, noise, SNRtag, SNRin_meas, SNRout, dSNR, RMSdec}
for i = 1:length(mixed_files)
    mixed_name = mixed_files(i).name;
    [~, base, ~] = fileparts(mixed_name);                

    parts = split(string(base), '_');                    
    if numel(parts) < 3, continue; end
    clean_base = char(lower(parts(1)));
    noise_tag  = char(parts(2));

    snr_tag = NaN;
    tok = regexp(base,'snr([+-]?\d+)dB','tokens','once');
    if ~isempty(tok), snr_tag = str2double(tok{1}); end
    snr_disp = 'NA'; if ~isnan(snr_tag), snr_disp = sprintf('%+d', snr_tag); end

    filtered_file = fullfile(filtered_path, [base '_wiener.wav']);
    if ~isfile(filtered_file) || ~isKey(cleanMap, clean_base), continue; end
    clean_file = cleanMap(clean_base);

    [x_clean, Fs_c] = audioread(clean_file);                      x_clean = mean(x_clean,2);
    [x_noisy, Fs_m] = audioread(fullfile(mixed_files(i).folder, mixed_name)); x_noisy = mean(x_noisy,2);
    [x_filt,  Fs_f] = audioread(filtered_file);                   x_filt  = mean(x_filt,2);
    if Fs_c ~= Fs_target, x_clean = resample(x_clean, Fs_target, Fs_c); end
    if Fs_m ~= Fs_target, x_noisy = resample(x_noisy, Fs_target, Fs_m); end
    if Fs_f ~= Fs_target, x_filt  = resample(x_filt,  Fs_target, Fs_f); end

    L = min([length(x_clean), length(x_noisy), length(x_filt)]);
    x_clean = x_clean(1:L); x_noisy = x_noisy(1:L); x_filt = x_filt(1:L);

    e_in  = x_noisy - x_clean;
    e_out = x_filt  - x_clean;
    snr_b = 10*log10(sum(x_clean.^2) / (sum(e_in.^2)  + eps));
    snr_a = 10*log10(sum(x_clean.^2) / (sum(e_out.^2) + eps));
    dSNR  = snr_a - snr_b;
    rms_b = sqrt(mean(e_in.^2));
    rms_a = sqrt(mean(e_out.^2));
    rms_reduce = (1 - rms_a/max(rms_b,eps))*100;

    fprintf('%-28s %-10s %-8s %-10.2f %-10.2f %-10.2f %-9.2f\n', ...
        mixed_name, noise_tag, snr_disp, snr_b, snr_a, dSNR, rms_reduce);

    M(end+1,:) = {mixed_name, noise_tag, snr_tag, snr_b, snr_a, dSNR, rms_reduce};
end
fprintf('----------------------------------------------------------------------------------------------\n Đánh giá hoàn tất\n');

% Save CSV
if saveCSV && ~isempty(M)
    T = cell2table(M, 'VariableNames', ...
        {'file','noise','SNRin_tag_dB','SNRin_meas_dB','SNRout_dB','DeltaSNR_dB','RMS_Reduction_pct'});
    if ~exist('outputs','dir'), mkdir('outputs'); end
    outCSV = fullfile('outputs','wiener_compare_metrics.csv');
    writetable(T, outCSV);
    fprintf('Đã lưu CSV: %s\n', outCSV);
end

% Nghe & so sánh 1 mẫu (DeltaSNR lớn nhất)
if (doFigure || doListen) && exist('T','var') && ~isempty(T)
    [~, idxBest] = max(T.DeltaSNR_dB);

    sampleFile = T.file{idxBest};             
    [~, base, ~] = fileparts(sampleFile);
    parts_base = split(string(base), '_');
    clean_key  = char(lower(parts_base(1)));

    f_clean = cleanMap(clean_key);
    f_noisy = fullfile(mixed_path, [base '.wav']);
    f_filt  = fullfile(filtered_path, [base '_wiener.wav']);

    if isfile(f_clean) && isfile(f_noisy) && isfile(f_filt)
        [xc, Fc] = audioread(f_clean);  xc = mean(xc,2);
        [xn, Fn] = audioread(f_noisy);  xn = mean(xn,2);
        [xf, Ff] = audioread(f_filt);   xf = mean(xf,2);

        if Fc ~= Fs_target, xc = resample(xc, Fs_target, Fc); end
        if Fn ~= Fs_target, xn = resample(xn, Fs_target, Fn); end
        if Ff ~= Fs_target, xf = resample(xf, Fs_target, Ff); end

        Ls = min([length(xc), length(xn), length(xf)]);
        xc = xc(1:Ls); xn = xn(1:Ls); xf = xf(1:Ls);

        if doListen
            try
                sound([xn; zeros(round(0.3*Fs_target),1); xf; zeros(round(0.3*Fs_target),1); xc], Fs_target);
            catch
                warning('Không phát được âm thanh trên hệ thống này.');
            end
        end

        if doFigure
            figure('Name',['So sánh: ' strrep(base,'_','\_')], 'NumberTitle','off','Position',[100 100 1000 650]);
            subplot(3,2,1); plot(xn); title('Noisy');   xlim([1 Ls]);
            subplot(3,2,3); plot(xf); title('Filtered');xlim([1 Ls]);
            subplot(3,2,5); plot(xc); title('Clean');   xlim([1 Ls]);
            subplot(3,2,2); spectrogram(xn,256,128,256,Fs_target,'yaxis'); title('Noisy');
            subplot(3,2,4); spectrogram(xf,256,128,256,Fs_target,'yaxis'); title('Filtered');
            subplot(3,2,6); spectrogram(xc,256,128,256,Fs_target,'yaxis'); title('Clean');
            sgtitle(['So sánh Wiener: ' strrep(base,'_','\_')]);
        end
    end
end

