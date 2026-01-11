function avg = process_results(raw, data, config)
    % PROCESS_RESULTS - Oblicza średnie i zapisuje wyniki (z nazwą pliku)
    
    avg.times_enc = mean(raw.times_enc, 2, 'omitnan');
    avg.times_dec = mean(raw.times_dec, 2, 'omitnan');
    avg.npcr = mean(raw.npcr, 2, 'omitnan');
    avg.uaci = mean(raw.uaci, 2, 'omitnan');
    avg.entropy = mean(raw.entropy, 2, 'omitnan');
    avg.corr_h = mean(raw.corr_h, 2, 'omitnan');
    avg.corr_v = mean(raw.corr_v, 2, 'omitnan');
    avg.corr_d = mean(raw.corr_d, 2, 'omitnan');

    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    
    % Wyciągamy samą nazwę pliku (bez ścieżki i rozszerzenia) do nazwy pliku wynikowego
    [~, fname_only, ~] = fileparts(data.filename);
    
    if ~exist('results', 'dir'), mkdir('results'); end
    
    % --- 1. Zapisz plik .MAT (z nazwą obrazu w nazwie pliku) ---
    filename_mat = sprintf('results/results_%s_%s.mat', fname_only, timestamp);
    
    titles = config.titles;
    N_RUNS = config.N_RUNS;
    
    save(filename_mat, 'titles', 'N_RUNS', 'raw', 'avg', 'data');
    fprintf('Wyniki (.MAT) zapisano w: %s\n', filename_mat);

    % --- 2. Zapisz plik .CSV (z nazwą obrazu w nazwie i kolumnie) ---
    filename_csv = sprintf('results/results_%s_%s.csv', fname_only, timestamp);
    
    fid = fopen(filename_csv, 'w');
    if fid == -1
        warning('Nie można utworzyć pliku CSV: %s', filename_csv);
        return;
    end
    
    % Nagłówek CSV - dodano kolumnę Image_File
    fprintf(fid, 'Image_File,Algorithm,Run_Index,Time_Enc_s,Time_Dec_s,NPCR_percent,UACI_percent,Entropy,Corr_Hor,Corr_Ver,Corr_Diag\n');
    
    [num_algos, num_runs] = size(raw.times_enc);
    
    for i = 1:num_algos
        algo_name = titles{i};
        % Usuwamy przecinki z nazwy algorytmu dla bezpieczeństwa CSV
        safe_name = strrep(algo_name, ',', ';'); 
        
        for r = 1:num_runs
            % Dane
            t_e = raw.times_enc(i, r);
            t_d = raw.times_dec(i, r);
            n   = raw.npcr(i, r);
            u   = raw.uaci(i, r);
            e   = raw.entropy(i, r);
            ch  = raw.corr_h(i, r);
            cv  = raw.corr_v(i, r);
            cd  = raw.corr_d(i, r);
            
            % Zapis wiersza z nazwą pliku na początku
            fprintf(fid, '%s,%s,%d,%.8f,%.8f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
                data.filename, safe_name, r, t_e, t_d, n, u, e, ch, cv, cd);
        end
    end
    
    fclose(fid);
    fprintf('Szczegóły (.CSV) zapisano w: %s\n', filename_csv);
end