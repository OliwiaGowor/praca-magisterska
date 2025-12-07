function visualize_results(avg, raw, data, config)
    % VISUALIZE_RESULTS - Table with normal notation (0.00...) and wide layout
    
    titles = config.titles;
    num_algos = length(titles);
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    
    % Prepare folder
    scriptPath = fileparts(mfilename('fullpath'));
    resultsFolder = fullfile(scriptPath, 'results');
    if ~exist(resultsFolder, 'dir'), mkdir(resultsFolder); end
    
    disp('Generowanie wizualizacji (Normalna notacja, szeroka tabela)...');

    %% ============================================================
    %% FIGURE 1: METRICS REPORT (Chart + Wide Table)
    %% ============================================================
    % Window width: 1500px (very wide)
    hFigMetrics = figure('Name', 'Raport Metryk', 'NumberTitle', 'off', ...
                         'Position', [50, 50, 1500, 850], 'Color', 'w');
    
    % --- 1. Top: Bar Chart ---
    subplot(2, 1, 1);
    set(gca, 'Position', [0.08, 0.50, 0.90, 0.45]); % Margins adjusted to width
    
    bar_data = [avg.times_enc avg.times_dec];
    b = bar(categorical(titles), bar_data);
    b(1).FaceColor = [0.2 0.6 0.8]; 
    b(2).FaceColor = [0.8 0.4 0.2]; 
    
    legend('Szyfrowanie', 'Deszyfrowanie', 'Location', 'best');
    title(sprintf('Porównanie Czasu Wykonania (Średnia z %d przebiegów)', config.N_RUNS), 'FontSize', 14);
    ylabel('Czas (sekundy)');
    grid on;
    
    % Logarithmic scale only if differences are extreme
    if max(avg.times_enc) / min(avg.times_enc(avg.times_enc>0)) > 50
        set(gca, 'YScale', 'log'); 
        ylabel('Czas (sekundy) - Skala Log');
    end

    % --- 2. Bottom: Data Table (Notation 0.000...) ---
    colNames = {'Algorytm', 'T_Enc [s]', 'T_Dec [s]', 'NPCR [%]', 'UACI [%]', 'Entropia', 'Korelacja'};
    tableData = cell(num_algos, 7);
    
    for i = 1:num_algos
        avg_corr = mean([avg.corr_h(i), avg.corr_v(i), avg.corr_d(i)], 'omitnan');
        
        % FORMATTING CHANGE:
        % %.8f -> Normal decimal notation (e.g., 0.00004505), 8 decimal places
        % %.5f -> Percentage metrics, 5 decimal places
        
        tableData{i,1} = titles{i};
        tableData{i,2} = sprintf('%.8f', avg.times_enc(i));  
        tableData{i,3} = sprintf('%.8f', avg.times_dec(i));
        tableData{i,4} = sprintf('%.5f', avg.npcr(i));       
        tableData{i,5} = sprintf('%.5f', avg.uaci(i));
        tableData{i,6} = sprintf('%.5f', avg.entropy(i));
        tableData{i,7} = sprintf('%.5f', avg_corr);
    end
    
    % Table height adapted to number of rows
    tableHeight = min(0.4, (num_algos + 2) * 0.045); 
    
    % Table stretched to 98% of window width
    t = uitable('Parent', hFigMetrics, 'Data', tableData, ...
                'ColumnName', colNames, ...
                'RowName', [], ...
                'Units', 'normalized', ...
                'Position', [0.01, 0.02, 0.98, tableHeight], ... 
                'FontSize', 11);
            
    % COLUMN WIDTHS (Sum ~1450px)
    % Name: 350px, Rest: 180px
    t.ColumnWidth = {350, 180, 180, 180, 180, 180, 180};

    % Save Figure 1
    file_metrics = fullfile(resultsFolder, sprintf('report_metrics_%s.png', timestamp));
    saveas(hFigMetrics, file_metrics);

    %% ============================================================
    %% FIGURE 2: VISUAL GALLERY
    %% ============================================================
    cols = 4; 
    rows = ceil(num_algos / 2) + 1; 
    
    hFigVisuals = figure('Name', 'Galeria Szyfrogramów', 'NumberTitle', 'off', ...
                         'Position', [150, 150, 1400, 250 * rows], 'Color', 'w');

    % Originals
    subplot(rows, cols, [2 3]); 
    imshow([data.img_orig, 255*ones(size(data.img_orig,1), 10, 'uint8'), data.img_mod]);
    title('LEWA: Oryginał  |  PRAWA: Zmodyfikowany (1px)');
    axis off;

    % Algorithms
    for i = 1:num_algos
        row_idx = ceil(i/2) + 1;
        is_second_in_row = mod(i, 2) == 0;
        if ~is_second_in_row, col_start = 1; else, col_start = 3; end
        
        if isempty(raw.images_enc{i}), continue; end
        
        subplot(rows, cols, (row_idx-1)*cols + col_start);
        imshow(raw.images_enc{i});
        title(sprintf('%s\n(Enc)', titles{i}), 'FontWeight', 'bold', 'Interpreter', 'none');
        
        subplot(rows, cols, (row_idx-1)*cols + col_start + 1);
        imshow(raw.images_dec{i});
        title('(Dec)');
    end

    % Save Figure 2
    file_visuals = fullfile(resultsFolder, sprintf('report_visuals_%s.png', timestamp));
    saveas(hFigVisuals, file_visuals);

    %% ============================================================
    %% CLEANUP
    %% ============================================================
    fprintf('--------------------------------------------------\n');
    fprintf(' Raport Metryk:   %s\n', file_metrics);
    fprintf(' Raport Wizualny: %s\n', file_visuals);
    fprintf('--------------------------------------------------\n');
    
    disp('Wizualizacja gotowa.');
end
