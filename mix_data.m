clear; clc;

% Options
rootDir   = 'data';
cleanDir  = fullfile(rootDir, 'clean');
noiseDir  = fullfile(rootDir, 'noise');
outDir    = fullfile(rootDir, 'mixed');
targetFs  = 16000;                % 16 kHz
snrLevels = [-5 0 5 10 15];       % dB
rng(0);                           

if ~exist(outDir, 'dir'); mkdir(outDir); end

cleanList = dir(fullfile(cleanDir, '*.flac'));
noiseList = dir(fullfile(noiseDir, '*.wav'));

assert(~isempty(cleanList), 'Không tìm thấy .flac trong data/clean/');
assert(~isempty(noiseList), 'Không tìm thấy .wav trong data/noise/');

nMade = 0;

for ic = 1:numel(cleanList)
    % Đọc clean (.flac)
    cleanPath = fullfile(cleanList(ic).folder, cleanList(ic).name);
    [x, fsX] = audioread(cleanPath);
    x = mean(x, 2);                           
    if fsX ~= targetFs
        x = resample(x, targetFs, fsX);
    end
    fs = targetFs;

    x = x ./ max(1e-12, local_rms(x));

    [~, cleanBase, ~] = fileparts(cleanList(ic).name);

    for in = 1:numel(noiseList)
        % Đọc noise (.wav)
        noisePath = fullfile(noiseList(in).folder, noiseList(in).name);
        [v, fsV] = audioread(noisePath);
        v = mean(v, 2);                      
        if fsV ~= fs
            v = resample(v, fs, fsV);
        end

        Lx = numel(x);
        if numel(v) < Lx
            rep = ceil(Lx/numel(v));
            v = repmat(v, rep, 1);
        end
        startIdx = randi(numel(v) - Lx + 1);
        vSegBase = v(startIdx:startIdx+Lx-1);

        [~, noiseBase, ~] = fileparts(noiseList(in).name);

        for is = 1:numel(snrLevels)
            SNRdB = snrLevels(is);

            vSeg = vSegBase;
            alpha = 10^(-SNRdB/20) * (local_rms(x) / max(1e-12, local_rms(vSeg)));
            vScaled = alpha * vSeg;

            y = x + vScaled;

            peak = max(abs(y));
            if peak > 0.999
                y = 0.999 * y / peak;
            end

            tag = sprintf('snr%+03ddB', round(SNRdB));
            outName = sprintf('%s_%s_%s.wav', cleanBase, noiseBase, tag);
            outPath = fullfile(outDir, outName);

            audiowrite(outPath, y, fs);
            nMade = nMade + 1;
            fprintf('[%04d] %s\n', nMade, outName);
        end
    end
end

fprintf('Hoàn tất: tạo %d file trong %s\n', nMade, outDir);

function r = local_rms(x)
    r = sqrt(mean(x.^2) + eps);
end
