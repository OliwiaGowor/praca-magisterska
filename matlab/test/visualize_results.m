function visualize_results(avg, raw, data, config, output_folder)
    % VISUALIZE_RESULTS - Generuje tabele PNG w podanym folderze
    
    titles = config.titles;
    num_algos = length(titles);
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    
    [~, fname_only, ~] = fileparts(data.filename);
    
    % Upewnij się, że folder istnieje
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end
    
    disp(['Generowanie wizualizacji w: ' output_folder]);

    %% --- TABELA 1: GŁÓWNA ---
    fig_height = max(400, 35*(num_algos+1) + 100);
    f1 = figure('Name', 'Glowne Metryki', 'Position', [50, 50, 1600, fig_height], 'Color', 'w', 'Visible', 'on', 'MenuBar', 'none');
    
    uicontrol('Style', 'text', 'String', ['Główne Metryki: ' fname_only], ...
              'Position', [20, fig_height-40, 1000, 25], 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
          
    colNames = {'Algorytm', 'Czas Szyfr. [s]', 'Czas Deszyfr. [s]', 'Entropia', 'Kor. Poz.', 'Kor. Pion.', 'Kor. Diag.', 'Kor. Śr.'};
    dat1 = cell(num_algos + 1, 8);
    
    % Oryginał
    if isfield(data, 'stats')
        e_orig = data.stats.entropy;
        c_h = data.stats.corr(1); c_v = data.stats.corr(2); c_d = data.stats.corr(3);
        c_avg = mean([c_h, c_v, c_d]);
        dat1(1,:) = {'>>> OBRAZ ORYGINALNY <<<', '-', '-', sprintf('%.5f', e_orig), sprintf('%.5f', c_h), sprintf('%.5f', c_v), sprintf('%.5f', c_d), sprintf('%.5f', c_avg)};
    else
        dat1(1,:) = {'Oryginał', '-', '-', '-', '-', '-', '-', '-'};
    end
    
    % Algorytmy
    for i=1:num_algos
        if i > length(avg.times_enc), continue; end
        c_avg_algo = mean([avg.corr_h(i), avg.corr_v(i), avg.corr_d(i)], 'omitnan');
        dat1(i+1,:) = {titles{i}, sprintf('%.6f', avg.times_enc(i)), sprintf('%.6f', avg.times_dec(i)), sprintf('%.5f', avg.entropy(i)), sprintf('%.5f', avg.corr_h(i)), sprintf('%.5f', avg.corr_v(i)), sprintf('%.5f', avg.corr_d(i)), sprintf('%.5f', c_avg_algo)};
    end
    
    uitable('Parent', f1, 'Data', dat1, 'ColumnName', colNames, 'Units', 'normalized', 'Position', [0.02, 0.02, 0.96, 0.88], 'ColumnWidth', {400, 120, 120, 100, 100, 100, 100, 100}, 'FontSize', 11);
    
    drawnow; frame1 = getframe(f1);
    imwrite(frame1.cdata, fullfile(output_folder, sprintf('table_main_%s.png', timestamp)));

    %% --- TABELA 2: WRAŻLIWOŚĆ ---
    f2 = figure('Name', 'Analiza Wrazliwosci', 'Position', [100, 100, 1800, max(400, 35*num_algos+100)], 'Color', 'w', 'Visible', 'on', 'MenuBar', 'none');
    uicontrol('Style', 'text', 'String', ['Analiza Wrażliwości (NPCR/UACI): ' fname_only], 'Position', [20, max(400, 35*num_algos+100)-40, 1000, 25], 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');

    colNames2 = {'Algorytm', 'NPCR (Pocz.)', 'UACI (Pocz.)', 'NPCR (Środek)', 'UACI (Środek)', 'NPCR (Koniec)', 'UACI (Koniec)'};
    dat2 = cell(num_algos, 7);
    for i=1:num_algos
        if i > length(avg.times_enc), continue; end
        dat2(i,:) = {titles{i}, sprintf('%.8f', avg.npcr_start(i)), sprintf('%.8f', avg.uaci_start(i)), sprintf('%.8f', avg.npcr_mid(i)), sprintf('%.8f', avg.uaci_mid(i)), sprintf('%.8f', avg.npcr_end(i)), sprintf('%.8f', avg.uaci_end(i))};
    end
    
    uitable('Parent', f2, 'Data', dat2, 'ColumnName', colNames2, 'Units', 'normalized', 'Position', [0.02, 0.02, 0.96, 0.88], 'FontSize', 10, 'ColumnWidth', {400, 130, 130, 130, 130, 130, 130});
    
    drawnow; frame2 = getframe(f2);
    imwrite(frame2.cdata, fullfile(output_folder, sprintf('table_sensitivity_%s.png', timestamp)));
    
    close(f1); close(f2);
end