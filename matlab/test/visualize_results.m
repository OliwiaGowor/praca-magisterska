function visualize_results(avg, raw, data, config)
    % VISUALIZE_RESULTS - Generuje wizualizacje z informacją o nazwie pliku
    
    titles = config.titles;
    num_algos = length(titles);
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    
    % Nazwa pliku do podpisów
    img_name = data.filename;
    [~, fname_only, ~] = fileparts(img_name);
    
    scriptPath = fileparts(mfilename('fullpath'));
    resultsFolder = fullfile(scriptPath, 'results');
    if ~exist(resultsFolder, 'dir'), mkdir(resultsFolder); end
    
    disp('Generowanie wizualizacji...');

    %% ============================================================
    %% 1. WYKRES (BAR CHART)
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
    
    % Tytuł z nazwą pliku
    title(sprintf('Czas Wykonania dla pliku: %s (Średnia z %d przebiegów)', img_name, config.N_RUNS), ...
          'FontSize', 16, 'Interpreter', 'none');
          
    ylabel('Czas (sekundy)', 'FontSize', 13);
    grid on;
    
    valid_times = t_enc(t_enc > 0);
    if ~isempty(valid_times) && (max(valid_times) / min(valid_times) > 50)
        set(gca, 'YScale', 'log'); 
        ylabel('Czas (sekundy) - Skala Log', 'FontSize', 13);
    end
    set(gca, 'Position', [0.06, 0.25, 0.93, 0.68]);

    % Zapis z nazwą pliku
    file_chart = fullfile(resultsFolder, sprintf('chart_%s_%s.png', fname_only, timestamp));
    saveas(hFigChart, file_chart);

    %% ============================================================
    %% 2. TABELA WYNIKÓW
    %% ============================================================
    row_height_px = 35; header_height_px = 40; font_size_table = 12; window_width = 1900;
    table_height_px = (num_algos * row_height_px) + header_height_px + 40; % +40 na tytuł
    
    hFigTable = figure('Name', 'Tabela Wyników', 'NumberTitle', 'off', ...
                       'Position', [50, 50, window_width, table_height_px], ...
                       'Color', 'w', 'Visible', 'on', 'MenuBar', 'none', 'ToolBar', 'none');

    % Tytuł tekstowy nad tabelą
    uicontrol('Style', 'text', 'String', sprintf('Wyniki dla pliku: %s', img_name), ...
              'Position', [20, table_height_px-30, 500, 25], ...
              'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', ...
              'BackgroundColor', 'w');

    colNames = {'Algorytm', 'T_Enc [s]', 'T_Dec [s]', 'NPCR [%]', 'UACI [%]', 'Entropia', ...
                'Corr H', 'Corr V', 'Corr D', 'Corr Avg'};
    
    tableData = cell(num_algos, 10);
    for i = 1:num_algos
        if i > length(avg.times_enc), continue; end
        avg_corr = mean([avg.corr_h(i), avg.corr_v(i), avg.corr_d(i)], 'omitnan');
        
        tableData{i,1} = titles{i};
        tableData{i,2} = sprintf('%.8f', avg.times_enc(i));  
        tableData{i,3} = sprintf('%.8f', avg.times_dec(i));
        tableData{i,4} = sprintf('%.5f', avg.npcr(i));       
        tableData{i,5} = sprintf('%.5f', avg.uaci(i));
        tableData{i,6} = sprintf('%.5f', avg.entropy(i));
        tableData{i,7} = sprintf('%.5f', avg.corr_h(i));
        tableData{i,8} = sprintf('%.5f', avg.corr_v(i));
        tableData{i,9} = sprintf('%.5f', avg.corr_d(i));
        tableData{i,10} = sprintf('%.5f', avg_corr);
    end
    
    t = uitable('Parent', hFigTable, 'Data', tableData, 'ColumnName', colNames, ...
                'RowName', [], 'Units', 'normalized', 'Position', [0, 0, 1, 0.94], ... 
                'FontSize', font_size_table);
    t.ColumnWidth = {450, 160, 160, 120, 120, 120, 120, 120, 120, 120};

    file_table = fullfile(resultsFolder, sprintf('table_%s_%s.png', fname_only, timestamp));
    drawnow; frame = getframe(hFigTable); imwrite(frame.cdata, file_table);

    %% ============================================================
    %% 3. GALERIA WIZUALNA
    %% ============================================================
    cols = 4; rows = ceil(num_algos / 2) + 1; 
    hFigVisuals = figure('Name', 'Galeria', 'NumberTitle', 'off', ...
                         'Position', [100, 100, 1600, 300 * rows], 'Color', 'w');

    subplot(rows, cols, [2 3]); 
    imshow([data.img_orig, 255*ones(size(data.img_orig,1), 10, 'uint8'), data.img_mod]);
    title(sprintf('Plik: %s (Oryginał vs Zmodyfikowany)', img_name), 'FontSize', 12, 'Interpreter', 'none');
    axis off;

    for i = 1:num_algos
        if i > length(raw.images_enc) || isempty(raw.images_enc{i}), continue; end
        row_idx = ceil(i/2) + 1;
        if mod(i, 2) ~= 0, col_start = 1; else, col_start = 3; end
        
        subplot(rows, cols, (row_idx-1)*cols + col_start);
        imshow(raw.images_enc{i});
        title(sprintf('%s\n(Enc)', titles{i}), 'FontWeight', 'bold', 'Interpreter', 'none', 'FontSize', 10);
        
        subplot(rows, cols, (row_idx-1)*cols + col_start + 1);
        imshow(raw.images_dec{i});
        title('(Dec)', 'FontSize', 10);
    end

    file_visuals = fullfile(resultsFolder, sprintf('visuals_%s_%s.png', fname_only, timestamp));
    saveas(hFigVisuals, file_visuals);

    % --- FIX: Odtworzenie nazwy pliku CSV dla raportu w konsoli ---
    % Zmienna filename_csv z process_results tu nie istnieje, więc tworzymy string ponownie
    csv_report_path = fullfile(resultsFolder, sprintf('results_%s_%s.csv', fname_only, timestamp));

    fprintf('--------------------------------------------------\n');
    fprintf(' Raporty dla pliku: %s\n', img_name);
    fprintf(' Wykres:  %s\n', file_chart);
    fprintf(' Tabela:  %s\n', file_table);
    fprintf(' Galeria: %s\n', file_visuals);
    fprintf(' CSV:     %s\n', csv_report_path);
    fprintf('--------------------------------------------------\n');
end