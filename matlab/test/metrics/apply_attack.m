function main_gui()
    % MAIN_GUI - Interfejs graficzny do testowania algorytmów
    
    % 1. KONFIGURACJA ŚCIEŻEK (Naprawia błąd "Undefined function")
    try
        setup_paths();
    catch
        addpath(genpath(fileparts(mfilename('fullpath')))); % Fallback
    end
    
    close all;
    
    % Lista algorytmów
    algo_list = { ...
        'AES-CBC (Hardware-C)', ...
        'Blowfish-CBC (Software-C)', ...
        'ChaCha20 (Software-C)', ...
        'Chaotyczny (MATLAB)', ...
        'Chaos Bit Shuffling (MATLAB)', ...
        'Hyperchaotic 2D (MATLAB)', ...
        'Hu & Tian - Two-Stage Logistic (MATLAB)', ...
        'Entropy Wavelet-Chaos (Pourasad et al.)'
    };

    % --- UI ---
    hFig = figure('Name', 'Analizator Szyfrowania', 'NumberTitle', 'off', ...
                  'Position', [100, 100, 1000, 700], 'MenuBar', 'none');

    % Panel Ładowania
    pnlLoad = uipanel('Parent', hFig, 'Title', 'Obraz', 'Position', [0.01, 0.55, 0.30, 0.44]);
    uicontrol('Parent', pnlLoad, 'Style', 'pushbutton', 'String', 'Załaduj Obraz', ...
              'Units', 'normalized', 'Position', [0.1, 0.85, 0.8, 0.1], 'Callback', @load_image_callback);
    axOrig = axes('Parent', pnlLoad, 'Position', [0.1, 0.1, 0.8, 0.7]);
    title(axOrig, 'Brak obrazu'); axis(axOrig, 'off');
    
    % Panel Benchmark
    pnlBench = uipanel('Parent', hFig, 'Title', 'Benchmark', 'Position', [0.01, 0.01, 0.30, 0.53]);
    btnBench = uicontrol('Parent', pnlBench, 'Style', 'pushbutton', 'String', 'URUCHOM TESTY', ...
              'Units', 'normalized', 'Position', [0.1, 0.4, 0.8, 0.15], 'Enable', 'off', ...
              'BackgroundColor', [0.8 1 0.8], 'Callback', @run_benchmark_callback);
    txtStatus = uicontrol('Parent', pnlBench, 'Style', 'text', 'String', 'Gotowy', ...
                          'Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.3]);

    % Panel Ataków
    pnlAttack = uipanel('Parent', hFig, 'Title', 'Testy Ataków', 'Position', [0.32, 0.01, 0.67, 0.98]);
    
    % Kontrolki
    uicontrol('Parent', pnlAttack, 'Style', 'text', 'String', 'Algorytm:', 'Units', 'normalized', 'Position', [0.02, 0.93, 0.15, 0.03]);
    popAlgo = uicontrol('Parent', pnlAttack, 'Style', 'popupmenu', 'String', algo_list, 'Units', 'normalized', 'Position', [0.18, 0.93, 0.30, 0.03]);
    
    uicontrol('Parent', pnlAttack, 'Style', 'text', 'String', 'Atak:', 'Units', 'normalized', 'Position', [0.50, 0.93, 0.10, 0.03]);
    popType = uicontrol('Parent', pnlAttack, 'Style', 'popupmenu', 'String', {'Crop', 'Salt & Pepper', 'Gaussian', 'Tamper'}, ...
                        'Units', 'normalized', 'Position', [0.60, 0.93, 0.15, 0.03]);
    
    uicontrol('Parent', pnlAttack, 'Style', 'text', 'String', 'Siła (0-1):', 'Units', 'normalized', 'Position', [0.02, 0.88, 0.20, 0.03]);
    edtInt = uicontrol('Parent', pnlAttack, 'Style', 'edit', 'String', '0.1', 'Units', 'normalized', 'Position', [0.22, 0.88, 0.10, 0.03]);
    
    btnAttack = uicontrol('Parent', pnlAttack, 'Style', 'pushbutton', 'String', 'Wykonaj Atak', ...
                          'Units', 'normalized', 'Position', [0.40, 0.87, 0.20, 0.05], 'Enable', 'off', 'Callback', @run_attack_callback);

    % Wykresy
    axEnc = axes('Parent', pnlAttack, 'Position', [0.05, 0.55, 0.40, 0.25]); title(axEnc, 'Zaszyfrowany'); axis(axEnc, 'off');
    axAtt = axes('Parent', pnlAttack, 'Position', [0.55, 0.55, 0.40, 0.25]); title(axAtt, 'Po Ataku'); axis(axAtt, 'off');
    axDec = axes('Parent', pnlAttack, 'Position', [0.30, 0.15, 0.40, 0.35]); title(axDec, 'Odszyfrowany'); axis(axDec, 'off');
    lblMetrics = uicontrol('Parent', pnlAttack, 'Style', 'text', 'String', '', 'Units', 'normalized', 'Position', [0.10, 0.05, 0.80, 0.05], 'FontSize', 12, 'FontWeight', 'bold');

    current_img = [];

    % --- CALLBACKS ---
    function load_image_callback(~, ~)
        [file, path] = uigetfile({'*.jpg;*.png;*.tif;*.bmp'}, 'Wybierz obraz');
        if isequal(file, 0), return; end
        try
            img = imread(fullfile(path, file));
            if size(img, 3) == 3, img = rgb2gray(img); end
            current_img = img;
            imshow(current_img, 'Parent', axOrig); title(axOrig, 'Oryginał');
            set(btnBench, 'Enable', 'on'); set(btnAttack, 'Enable', 'on');
            set(txtStatus, 'String', 'Obraz wczytany.');
        catch ME
            errordlg(ME.message, 'Błąd pliku');
        end
    end

    function run_benchmark_callback(~, ~)
        if isempty(current_img), return; end
        set(hFig, 'Pointer', 'watch'); set(txtStatus, 'String', 'Benchmark w toku...'); drawnow;
        try
            % --- Test configuration ---
            config.N_RUNS = 2; 

            % List of algorithm titles
            config.titles = { ...
                'AES-CBC (Hardware-C)', ...
                'Blowfish-CBC (Software-C)', ...
                'ChaCha20 (Software-C)', ...
                'Chaotyczny (MATLAB)', ...
                'DES-CBC (MATLAB)', ...
                'AES-CBC (Pure MATLAB)', ...
                'Blowfish-CBC (Pure MATLAB)', ...
                'ChaCha20 (Pure MATLAB)', ...
                '3D-Logistic-Chirikov (MATLAB)', ...
                'Chaos Bit Shuffling (MATLAB)', ...
                'Hyperchaotic 2D (MATLAB*)', ...
                'Hu & Tian - Two-Stage Logistic (MATLAB*)', ...
                };

            

            % 1. Data preparation
            [data] = setup_test_data(current_img);

            % 2. Run benchmarks
            raw_results = run_benchmarks(data, config);

            % 3. Process and save results
            avg_results = process_results(raw_results, data, config);

            visualize_results(avg_results, raw_results, data, config);
            set(txtStatus, 'String', 'Zakończono.');
        catch ME
            errordlg(ME.message, 'Błąd Benchmarku');
            set(txtStatus, 'String', 'Błąd.');
            disp(getReport(ME)); % Wypisz szczegóły w konsoli
        end
        set(hFig, 'Pointer', 'arrow');
    end

    function run_attack_callback(~, ~)
        if isempty(current_img), return; end
        idx = get(popAlgo, 'Value'); name = algo_list{idx};
        types = {'crop', 'salt_pepper', 'gaussian', 'tamper'};
        atype = types{get(popType, 'Value')};
        val = str2double(get(edtInt, 'String'));
        
        set(hFig, 'Pointer', 'watch'); set(lblMetrics, 'String', 'Przetwarzanie...'); drawnow;
        try
            [enc, dec] = get_algorithm_handles(name);
            
            % Szyfrowanie
            [C, k] = enc(current_img);
            imshow(uint8(C), 'Parent', axEnc); title(axEnc, 'Zaszyfrowany');
            
            % Atak
            tamper_img = [];
            if strcmp(atype, 'tamper')
                 tamper_img = uint8(randi([0 255], size(current_img))); % Szum jako intruz
            end
            C_att = apply_attack(uint8(C), atype, val, tamper_img);
            imshow(C_att, 'Parent', axAtt); title(axAtt, ['Atak: ' atype]);
            
            % Deszyfrowanie
            P_rec = dec(C_att, k);
            if ~isequal(size(P_rec), size(current_img)), P_rec = imresize(P_rec, size(current_img)); end
            imshow(uint8(P_rec), 'Parent', axDec); title(axDec, 'Odzyskany');
            
            % Metryki
            p = psnr(uint8(P_rec), current_img);
            s = ssim(uint8(P_rec), current_img);
            set(lblMetrics, 'String', sprintf('PSNR: %.2f dB, SSIM: %.4f', p, s));
            
        catch ME
            errordlg(ME.message, 'Błąd Ataku');
            set(lblMetrics, 'String', 'Błąd');
        end
        set(hFig, 'Pointer', 'arrow');
    end
end