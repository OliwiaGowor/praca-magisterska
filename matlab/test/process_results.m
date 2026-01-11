function avg = process_results(raw, data, config)
    % PROCESS_RESULTS - Zapisuje wyniki (w tym ORYGINAŁ w CSV)
    
    % Średnie
    avg.times_enc = mean(raw.times_enc, 2, 'omitnan');
    avg.times_dec = mean(raw.times_dec, 2, 'omitnan');
    avg.entropy   = mean(raw.entropy, 2, 'omitnan');
    avg.corr_h    = mean(raw.corr_h, 2, 'omitnan');
    avg.corr_v    = mean(raw.corr_v, 2, 'omitnan');
    avg.corr_d    = mean(raw.corr_d, 2, 'omitnan');
    
    avg.npcr_start = mean(raw.npcr_start, 2, 'omitnan'); avg.uaci_start = mean(raw.uaci_start, 2, 'omitnan');
    avg.npcr_mid   = mean(raw.npcr_mid, 2, 'omitnan');   avg.uaci_mid   = mean(raw.uaci_mid, 2, 'omitnan');
    avg.npcr_end   = mean(raw.npcr_end, 2, 'omitnan');   avg.uaci_end   = mean(raw.uaci_end, 2, 'omitnan');

    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    [~, fname_only, ~] = fileparts(data.filename);
    if ~exist('results', 'dir'), mkdir('results'); end
    
    save(sprintf('results/full_results_%s_%s.mat', fname_only, timestamp), 'avg', 'raw', 'config');

    % --- CSV GŁÓWNY ---
    f_main = sprintf('results/main_metrics_%s_%s.csv', fname_only, timestamp);
    fid = fopen(f_main, 'w');
    fprintf(fid, 'Plik,Algorytm,Przebieg,Czas_Szyfr,Czas_Deszyfr,Entropia,Kor_Poz,Kor_Pion,Kor_Diag,Kor_Srednia\n');
    
    % 1. Wiersz z ORYGINAŁEM
    if isfield(data, 'stats')
        e_orig = data.stats.entropy;
        c_h = data.stats.corr(1); c_v = data.stats.corr(2); c_d = data.stats.corr(3);
        c_avg = mean([c_h, c_v, c_d]);
        fprintf(fid, '%s,ORIGINAL_IMAGE,0,0,0,%.8f,%.8f,%.8f,%.8f,%.8f\n', ...
            data.filename, e_orig, c_h, c_v, c_d, c_avg);
    end
    
    % 2. Wiersze z ALGORYTMAMI
    [num_algos, num_runs] = size(raw.times_enc);
    for i=1:num_algos
        safe_name = strrep(config.titles{i}, ',', ';');
        for r=1:num_runs
            c_avg_algo = mean([raw.corr_h(i,r), raw.corr_v(i,r), raw.corr_d(i,r)]);
            fprintf(fid, '%s,%s,%d,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n', ...
                data.filename, safe_name, r, ...
                raw.times_enc(i,r), raw.times_dec(i,r), raw.entropy(i,r), ...
                raw.corr_h(i,r), raw.corr_v(i,r), raw.corr_d(i,r), c_avg_algo);
        end
    end
    fclose(fid);
    
    % --- CSV WRAŻLIWOŚĆ ---
    f_sens = sprintf('results/npcr_uaci_analysis_%s_%s.csv', fname_only, timestamp);
    fid = fopen(f_sens, 'w');
    fprintf(fid, 'Plik,Algorytm,Przebieg,NPCR_Poczatek,UACI_Poczatek,NPCR_Srodek,UACI_Srodek,NPCR_Koniec,UACI_Koniec\n');
    
    for i=1:num_algos
        safe_name = strrep(config.titles{i}, ',', ';');
        for r=1:num_runs
            fprintf(fid, '%s,%s,%d,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n', ...
                data.filename, safe_name, r, ...
                raw.npcr_start(i,r), raw.uaci_start(i,r), ...
                raw.npcr_mid(i,r),   raw.uaci_mid(i,r), ...
                raw.npcr_end(i,r),   raw.uaci_end(i,r));
        end
    end
    fclose(fid);
    
    fprintf('CSV Główny (z oryginałem): %s\n', f_main);
end