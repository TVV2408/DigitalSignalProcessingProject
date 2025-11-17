clear; clc; close all;

% Options
Fs_target   = 16000;
rootClean   = 'data/clean';
rootMixed   = 'data/mixed';
rootWiener  = 'outputs/wiener_filtered_fix';
rootLMS     = 'outputs/lms_filtered_fix';

% Chọn mẫu
speech_id = 'speech07';   % 'speech04'
noise_tag = 'ac';         % 'ac' | 'fan' | 'office'
snr_db    = 0;            % -5, 0, 5, 10, 15

if snr_db >= 0
    snr_str = sprintf('snr+%02ddB', snr_db);
else
    snr_str = sprintf('snr-%02ddB', abs(snr_db));
end
base = sprintf('%s_%s_%s', speech_id, noise_tag, snr_str);

f_clean  = fullfile(rootClean,  [speech_id '.flac']);
f_mixed  = fullfile(rootMixed,  [base '.wav']);
f_wiener = fullfile(rootWiener, [base '_wiener.wav']);
f_lms    = fullfile(rootLMS,    [base '_lms.wav']);

assert(isfile(f_clean),  'Không thấy clean: %s',  f_clean);
assert(isfile(f_mixed),  'Không thấy mixed: %s',  f_mixed);
assert(isfile(f_wiener), 'Không thấy wiener: %s', f_wiener);
assert(isfile(f_lms),    'Không thấy lms: %s',    f_lms);

[xc, Fs_c] = audioread(f_clean);   xc = mean(xc,2);
[xn, Fs_n] = audioread(f_mixed);   xn = mean(xn,2);
[xw, Fs_w] = audioread(f_wiener);  xw = mean(xw,2);
[xl, Fs_l] = audioread(f_lms);     xl = mean(xl,2);

if Fs_c ~= Fs_target, xc = resample(xc, Fs_target, Fs_c); end
if Fs_n ~= Fs_target, xn = resample(xn, Fs_target, Fs_n); end
if Fs_w ~= Fs_target, xw = resample(xw, Fs_target, Fs_w); end
if Fs_l ~= Fs_target, xl = resample(xl, Fs_target, Fs_l); end
Fs = Fs_target;

L = min([length(xc), length(xn), length(xw), length(xl)]);
xc = xc(1:L);  xn = xn(1:L);  xw = xw(1:L);  xl = xl(1:L);

snr_b = 10*log10(sum(xc.^2) / (sum((xn - xc).^2) + eps));
snr_w = 10*log10(sum(xc.^2) / (sum((xw - xc).^2) + eps));
snr_l = 10*log10(sum(xc.^2) / (sum((xl - xc).^2) + eps));

rms_w = sqrt(mean((xw - xc).^2));
rms_l = sqrt(mean((xl - xc).^2));

fprintf('\n SO SÁNH WIENER vs LMS – %s\n', base);
fprintf('------------------------------------------------------\n');
fprintf('SNR trước lọc  : %6.2f dB\n', snr_b);
fprintf('SNR Wiener     : %6.2f dB (Δ = %+6.2f)\n', snr_w, snr_w - snr_b);
fprintf('SNR LMS        : %6.2f dB (Δ = %+6.2f)\n', snr_l, snr_l - snr_b);
fprintf('RMS Wiener     : %.6f\n', rms_w);
fprintf('RMS LMS        : %.6f\n', rms_l);
fprintf('------------------------------------------------------\n');

try
    gap = zeros(round(0.3*Fs),1);
    sound([xn; gap; xw; gap; xl; gap; xc], Fs);
catch
    warning('Không phát được âm thanh trên hệ thống này.');
end

% ===== Vẽ biểu đồ =====
figure('Name', ['So sánh Wiener & LMS – ' strrep(base,'_','\_')], ...
       'NumberTitle','off','Position',[100 100 1000 650]);

subplot(2,2,[1 2]); hold on; grid on;
plot(xn, 'Color', [0.6 0.6 0.6], 'DisplayName', 'Noisy');
plot(xw, 'r',  'DisplayName', 'Wiener');
plot(xl, 'b',  'DisplayName', 'LMS');
plot(xc, 'k--','DisplayName', 'Clean');
xlabel('Mẫu (Sample)'); ylabel('Biên độ');
title('Waveform so sánh'); legend('Location','northeastoutside'); hold off;

% Spectrogram Wiener
subplot(2,2,3);
spectrogram(xw, 256, 128, 256, Fs, 'yaxis');
title('Wiener – Spectrogram');

% Spectrogram LMS
subplot(2,2,4);
spectrogram(xl, 256, 128, 256, Fs, 'yaxis');
title('LMS – Spectrogram');

sgtitle(['So sánh hai bộ lọc trên ' strrep(base,'_','\_')]);
