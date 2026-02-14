function visualize_results(avg, raw, data, config)
    % VISUALIZE_RESULTS - Generates a full set of visualizations:
    % 1. Time bar chart
    % 2. Main Table 
    % 3. Sensitivity Table (NPCR/UACI)
    % 4. Two visual galleries
    
    titles = config.titles;
    num_algos = length(titles);
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    
    if isfield(data, 'filename')
        [~, fname_only, ~] = fileparts(data.filename);
        full_filename = data.filename;
    else
        fname_only = 'unknown';
        full_filename = 'Nieznany Plik';
    end
    
    scriptPath = fileparts(mfilename('fullpath'));
    resultsFolder = fullfile(scriptPath, 'results');
    if ~exist(resultsFolder, 'dir'), mkdir(resultsFolder); end
    
    disp('Generowanie wizualizacji...');

    %% ============================================================
    %% 1. CHART (BAR CHART)
    %% ============================================================
    hFigChart = figure('Name', 'Wykres Czasu', 'NumberTitle', 'off', ...
                       'Position', [50, 50, 1600, 900], 'Color', 'w', 'Visible', 'on');
    
    t_enc = avg.times_enc(:);
    t_dec = avg.times_dec(:);
    min_len = min([length(t_enc), length(t_dec), length(titles)]);
    bar_data = [t_enc(1:min_len), t_dec(1:min_len)];
    
    b = bar(bar_data);
    b(1).FaceColor = [0.2 0.6 0.8]; 
    b(2).FaceColor = [0.8 0.4 0.2]; 
    
    set(gca, 'XTick', 1:min_len, 'XTickLabel', titles(1:min_len), 'FontSize', 10);
    set(gca, 'TickLabelInterpreter', 'none'); 
    xtickangle(45); 
    legend('Szyfrowanie', 'Deszyfrowanie', 'Location', 'northeast', 'FontSize', 12);
    
    title(sprintf('Czas Wykonania: %s (Średnia z %d przebiegów)', full_filename, config.N_RUNS), ...
          'FontSize', 16, 'Interpreter', 'none');
    ylabel('Czas (sekundy)', 'FontSize', 13);
    grid on;
    
    valid_times = t_enc(t_enc > 0);
    if ~isempty(valid_times) && (max(valid_times) / min(valid_times) > 50)
        set(gca, 'YScale', 'log'); 
        ylabel('Czas (sekundy) - Skala Log', 'FontSize', 13);
    end
    set(gca, 'Position', [0.06, 0.25, 0.93, 0.68]);

    file_chart = fullfile(resultsFolder, sprintf('chart_%s_%s.png', fname_only, timestamp));
    saveas(hFigChart, file_chart);

    %% ============================================================
    %% 2. TABLE 1: MAIN
    %% ============================================================
    fig_height = max(400, 35*(num_algos+1) + 100);
    f1 = figure('Name', 'Glowne Metryki', 'Position', [50, 50, 1600, fig_height], ...
                'Color', 'w', 'Visible', 'on', 'MenuBar', 'none');
    
    uicontrol('Style', 'text', 'String', ['Główne Metryki: ' full_filename], ...
              'Position', [20, fig_height-40, 1000, 25], ...
              'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
          
    colNames = {'Algorytm', 'Czas Szyfr. [s]', 'Czas Deszyfr. [s]', 'Entropia', 'Kor. Poz.', 'Kor. Pion.', 'Kor. Diag.', 'Kor. Śr.'};
    dat1 = cell(num_algos + 1, 8);
    
    % Original Row
    if isfield(data, 'stats')
        c_avg = mean([data.stats.corr(1), data.stats.corr(2), data.stats.corr(3)]);
        dat1(1,:) = {'>>> OBRAZ ORYGINALNY <<<', '-', '-', sprintf('%.5f', data.stats.entropy), ...
                     sprintf('%.5f', data.stats.corr(1)), sprintf('%.5f', data.stats.corr(2)), ...
                     sprintf('%.5f', data.stats.corr(3)), sprintf('%.5f', c_avg)};
    else
        dat1(1,:) = {'Oryginał', '-', '-', '-', '-', '-', '-', '-'};
    end
    
    for i=1:num_algos
        if i > length(avg.times_enc), continue; end
        c_avg_algo = mean([avg.corr_h(i), avg.corr_v(i), avg.corr_d(i)], 'omitnan');
        dat1(i+1,:) = {titles{i}, sprintf('%.6f', avg.times_enc(i)), sprintf('%.6f', avg.times_dec(i)), ...
                      sprintf('%.5f', avg.entropy(i)), sprintf('%.5f', avg.corr_h(i)), ...
                      sprintf('%.5f', avg.corr_v(i)), sprintf('%.5f', avg.corr_d(i)), sprintf('%.5f', c_avg_algo)};
    end
    
    uitable('Parent', f1, 'Data', dat1, 'ColumnName', colNames, ...
            'Units', 'normalized', 'Position', [0.02, 0.02, 0.96, 0.88], ...
            'ColumnWidth', {400, 120, 120, 100, 100, 100, 100, 100}, 'FontSize', 11);
    
    drawnow; frame1 = getframe(f1);
    file_table_main = fullfile(resultsFolder, sprintf('table_main_%s_%s.png', fname_only, timestamp));
    imwrite(frame1.cdata, file_table_main);

    %% ============================================================
    %% 3. TABLE 2: SENSITIVITY
    %% ============================================================
    f2 = figure('Name', 'Analiza Wrazliwosci', 'Position', [100, 100, 1800, max(400, 35*num_algos+100)], ...
                'Color', 'w', 'Visible', 'on', 'MenuBar', 'none');
    uicontrol('Style', 'text', 'String', ['Analiza Wrażliwości (NPCR/UACI): ' full_filename], ...
              'Position', [20, max(400, 35*num_algos+100)-40, 1000, 25], ...
              'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
    colNames2 = {'Algorytm', 'NPCR (Pocz.)', 'UACI (Pocz.)', 'NPCR (Środek)', 'UACI (Środek)', 'NPCR (Koniec)', 'UACI (Koniec)'};
    dat2 = cell(num_algos, 7);
    for i=1:num_algos
        if i > length(avg.times_enc), continue; end
        dat2(i,:) = {titles{i}, sprintf('%.8f', avg.npcr_start(i)), sprintf('%.8f', avg.uaci_start(i)), ...
            sprintf('%.8f', avg.npcr_mid(i)), sprintf('%.8f', avg.uaci_mid(i)), ...
            sprintf('%.8f', avg.npcr_end(i)), sprintf('%.8f', avg.uaci_end(i))};
    end
    uitable('Parent', f2, 'Data', dat2, 'ColumnName', colNames2, 'Units', 'normalized', ...
            'Position', [0.02, 0.02, 0.96, 0.88], 'FontSize', 10, 'ColumnWidth', {400, 130, 130, 130, 130, 130, 130});
    drawnow; frame2 = getframe(f2);
    file_table_sens = fullfile(resultsFolder, sprintf('table_sensitivity_%s_%s.png', fname_only, timestamp));
    imwrite(frame2.cdata, file_table_sens);

    %% ============================================================
    %% 4. VISUAL GALLERIES
    %% ============================================================
    
    % Algorithm categorization
    classic_idxs = [];
    chaos_idxs = [];
    
    for i = 1:num_algos
        t_lower = lower(titles{i});
        % If name contains classic keywords
        if contains(t_lower, 'des') || contains(t_lower, 'aes') || ...
           contains(t_lower, 'blowfish') || contains(t_lower, 'chacha')
            classic_idxs(end+1) = i;
        else
            chaos_idxs(end+1) = i;
        end
    end
    
    % Helper for gallery generation
    generate_gallery(classic_idxs, 'KLASYCZNE', 'visuals_classic', raw, titles, data, resultsFolder, fname_only, timestamp);
    generate_gallery(chaos_idxs, 'CHAOTYCZNE', 'visuals_chaos', raw, titles, data, resultsFolder, fname_only, timestamp);

    %% ============================================================
    %% SUMMARY
    %% ============================================================
    fprintf('--------------------------------------------------\n');
    fprintf(' Wygenerowano raporty dla pliku: %s\n', full_filename);
    fprintf(' Wykres:         %s\n', file_chart);
    fprintf(' Tabele:         %s, %s\n', file_table_main, file_table_sens);
    fprintf(' Galeria Klas.: %s\n', fullfile(resultsFolder, sprintf('visuals_classic_%s_%s.png', fname_only, timestamp)));
    fprintf(' Galeria Chaos: %s\n', fullfile(resultsFolder, sprintf('visuals_chaos_%s_%s.png', fname_only, timestamp)));
    fprintf('--------------------------------------------------\n');
end

function generate_gallery(indices, label, prefix, raw, titles, data, folder, fname, ts)
    if isempty(indices), return; end
    
    cols = 2; 

    rows_count = length(indices) + 1;
    
    fig_h = max(400, 300 * rows_count);
    
    hFig = figure('Name', ['Galeria ' label], 'NumberTitle', 'off', ...
                  'Position', [50, 50, 1200, fig_h], 'Color', 'w', 'Visible', 'on');
    
    % --- Row 1: Original and Modified ---
    subplot(rows_count, cols, 1);
    if isfield(data, 'img_orig'), imshow(data.img_orig); end
    title('Oryginał', 'FontSize', 12, 'FontWeight', 'bold');
    
    subplot(rows_count, cols, 2);
    if isfield(data, 'img_mod'), imshow(data.img_mod); end
    title('Zmodyfikowany (dla testów wrażliwości)', 'FontSize', 10);
    
    % --- Subsequent rows: Algorithms ---
    for k = 1:length(indices)
        algo_idx = indices(k);
        current_row = k + 1;
        
        % Column 1: Ciphertext
        subplot(rows_count, cols, (current_row-1)*cols + 1);
        if ~isempty(raw.images_enc{algo_idx})
            imshow(raw.images_enc{algo_idx});
        end
        title(sprintf('%s (Enc)', titles{algo_idx}), 'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');
        
        % Column 2: Decrypted
        subplot(rows_count, cols, (current_row-1)*cols + 2);
        if ~isempty(raw.images_dec{algo_idx})
            imshow(raw.images_dec{algo_idx});
        end
        title('(Dec)', 'FontSize', 10);
    end
    
    file_out = fullfile(folder, sprintf('%s_%s_%s.png', prefix, fname, ts));
    saveas(hFig, file_out);
end
