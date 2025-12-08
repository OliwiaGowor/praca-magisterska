function visualize_results(avg, raw, data, config)
    % VISUALIZE_RESULTS - Generuje wykresy i tabele (Fixed for uitable error)
    
    titles = config.titles;
    num_algos = length(titles);
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    
    % Prepare folder
    scriptPath = fileparts(mfilename('fullpath'));
    resultsFolder = fullfile(scriptPath, 'results');
    if ~exist(resultsFolder, 'dir'), mkdir(resultsFolder); end
    
    disp('Generowanie wizualizacji (Fix: uitable capture)...');

    %% ============================================================
    %% FIGURE 1: METRICS REPORT (Chart + Wide Table)
    %% ============================================================
    % Ustawiamy 'Visible', 'on' aby getframe zadziałało
    hFigMetrics = figure('Name', 'Raport Metryk', 'NumberTitle', 'off', ...
                         'Position', [50, 50, 1500, 850], 'Color', 'w', 'Visible', 'on');
    
    % --- 1. Top: Bar Chart ---
    subplot(2, 1, 1);
    set(gca, 'Position', [0.08, 0.50, 0.90, 0.45]); 
    
    t_enc = avg.times_enc(:);
    t_dec = avg.times_dec(:);
    
    % Zabezpieczenie wymiarów
    min_len = min([length(t_enc), length(t_dec), length(titles)]);
    if length(t_enc) ~= length(titles)
        warning('Liczba wyników (%d) różni się od liczby tytułów (%d). Przycinam.', length(t_enc), length(titles));
    end
    
    bar_data = [t_enc(1:min_len), t_dec(1:min_len)];
    plot_titles = titles(1:min_len);
    
    b = bar(bar_data);
    b(1).FaceColor = [0.2 0.6 0.8]; 
    b(2).FaceColor = [0.8 0.4 0.2]; 
    
    set(gca, 'XTick', 1:min_len, 'XTickLabel', plot_titles);
    xtickangle(45); 
    
    legend('Szyfrowanie', 'Deszyfrowanie', 'Location', 'best');
    title(sprintf('Porównanie Czasu Wykonania (Średnia z %d przebiegów)', config.N_RUNS), 'FontSize', 14);
    ylabel('Czas (sekundy)');
    grid on;
    
    valid_times = t_enc(t_enc > 0);
    if ~isempty(valid_times) && (max(valid_times) / min(valid_times) > 50)
        set(gca, 'YScale', 'log'); 
        ylabel('Czas (sekundy) - Skala Log');
    end

    % --- 2. Bottom: Data Table ---
    colNames = {'Algorytm', 'T_Enc [s]', 'T_Dec [s]', 'NPCR [%]', 'UACI [%]', 'Entropia', 'Korelacja'};
    tableData = cell(num_algos, 7);
    
    for i = 1:num_algos
        if i > length(avg.times_enc)
            tableData{i,1} = titles{i};
            tableData{i,2} = 'N/A';
            continue;
        end

        avg_corr = mean([avg.corr_h(i), avg.corr_v(i), avg.corr_d(i)], 'omitnan');
        
        tableData{i,1} = titles{i};
        tableData{i,2} = sprintf('%.8f', avg.times_enc(i));  
        tableData{i,3} = sprintf('%.8f', avg.times_dec(i));
        tableData{i,4} = sprintf('%.5f', avg.npcr(i));       
        tableData{i,5} = sprintf('%.5f', avg.uaci(i));
        tableData{i,6} = sprintf('%.5f', avg.entropy(i));
        tableData{i,7} = sprintf('%.5f', avg_corr);
    end
    
    tableHeight = min(0.4, (num_algos + 2) * 0.045); 
    t = uitable('Parent', hFigMetrics, 'Data', tableData, ...
                'ColumnName', colNames, ...
                'RowName', [], ...
                'Units', 'normalized', ...
                'Position', [0.01, 0.02, 0.98, tableHeight], ... 
                'FontSize', 11);
            
    t.ColumnWidth = {350, 180, 180, 180, 180, 180, 180};

    % --- FIX: Zapisywanie figury z tabelą (uitable) ---
    file_metrics = fullfile(resultsFolder, sprintf('report_metrics_%s.png', timestamp));
    
    % Wymuszamy odświeżenie grafiki
    drawnow; 
    
    try
        % Metoda 1: getframe (działa zawsze, jakość ekranowa)
        frame = getframe(hFigMetrics);
        imwrite(frame.cdata, file_metrics);
    catch
        % Fallback: próba użycia exportgraphics (nowsze MATLAB-y), jeśli getframe zawiedzie
        try
            exportgraphics(hFigMetrics, file_metrics);
        catch
            warning('Nie udało się zapisać raportu z tabelą. Pomińmy ten krok.');
        end
    end

    %% ============================================================
    %% FIGURE 2: VISUAL GALLERY
    %% ============================================================
    cols = 4; 
    rows = ceil(num_algos / 2) + 1; 
    
    hFigVisuals = figure('Name', 'Galeria Szyfrogramów', 'NumberTitle', 'off', ...
                         'Position', [150, 150, 1400, 250 * rows], 'Color', 'w');

    subplot(rows, cols, [2 3]); 
    imshow([data.img_orig, 255*ones(size(data.img_orig,1), 10, 'uint8'), data.img_mod]);
    title('LEWA: Oryginał  |  PRAWA: Zmodyfikowany (1px)');
    axis off;

    for i = 1:num_algos
        if i > length(raw.images_enc) || isempty(raw.images_enc{i})
            continue; 
        end
        
        row_idx = ceil(i/2) + 1;
        is_second_in_row = mod(i, 2) == 0;
        if ~is_second_in_row, col_start = 1; else, col_start = 3; end
        
        subplot(rows, cols, (row_idx-1)*cols + col_start);
        imshow(raw.images_enc{i});
        title(sprintf('%s\n(Enc)', titles{i}), 'FontWeight', 'bold', 'Interpreter', 'none');
        
        subplot(rows, cols, (row_idx-1)*cols + col_start + 1);
        imshow(raw.images_dec{i});
        title('(Dec)');
    end

    file_visuals = fullfile(resultsFolder, sprintf('report_visuals_%s.png', timestamp));
    saveas(hFigVisuals, file_visuals); % Tutaj saveas jest bezpieczne (brak uitable)

    fprintf('--------------------------------------------------\n');
    fprintf(' Raport Metryk:   %s\n', file_metrics);
    fprintf(' Raport Wizualny: %s\n', file_visuals);
    fprintf('--------------------------------------------------\n');
    disp('Wizualizacja gotowa.');
end
