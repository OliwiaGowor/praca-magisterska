function data = setup_test_data()
    % SETUP_TEST_DATA - Wczytuje obraz, konwertuje na szary i przygotowuje dane
    
    % --- KONFIGURACJA PLIKU ---
    % target_file = 'lena.png'; 
    target_file = 'images/lena_gray_256.tif'; % Twoje duże zdjęcie
    
    % Sprawdzenie czy plik istnieje
    if exist(target_file, 'file')
        filename = target_file;
        img_orig = imread(filename);
    else
        warning('Plik %s nie istnieje. Używam "cameraman.tif" jako fallback.', target_file);
        filename = 'cameraman.tif';
        img_orig = imread(filename);
    end
    
    % --- FIX: KONWERSJA NA SKALĘ SZAROŚCI ---
    % Jeśli obraz jest kolorowy (3 kanały RGB), zamień go na 2D (Grayscale).
    % To rozwiązuje problemy z wymiarami przy łączeniu macierzy i upraszcza szyfrowanie.
    if size(img_orig, 3) == 3
        img_orig = rgb2gray(img_orig);
        fprintf('Info: Obraz "%s" został przekonwertowany na skalę szarości.\n', filename);
    end
    
    % --- Modified image (for NPCR/UACI) ---
    img_mod = img_orig;
    img_mod(1) = bitxor(img_mod(1), 1); % Change 1 bit
    
    
    % Analysis of original
    entropy_orig = calculate_entropy(img_orig);
    [ch, cv, cd] = calculate_correlation(img_orig);
    
    % Pack into output structure
    data.img_orig = img_orig;
    data.stats.entropy = entropy_orig;
    data.stats.corr = [ch, cv, cd];
    data.img_mod  = img_mod;
    data.size     = size(img_orig);
    data.filename = filename;
    
    % Spłaszczone wersje (teraz wektor będzie 1x(WxH), bez mnożnika x3)
    data.img_flat = img_orig(:)';      
    data.img_flat_mod = img_mod(:)';   
end
