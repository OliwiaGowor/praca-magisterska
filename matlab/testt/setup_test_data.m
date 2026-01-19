function data = setup_test_data()
    % SETUP_TEST_DATA - Przygotowuje obrazy (Oryginał + 3 Modyfikacje) + Statystyki
    
    % --- Konfiguracja ścieżki do pliku ---
    target_file = 'images/lena_gray_256.tif'; 
    % Możesz tu zmienić na: 'images/pl23_con00216.jpg' itp.
    
    if exist(target_file, 'file')
        filename = target_file;
        img_orig = imread(filename);
    else
        warning('Plik %s nie istnieje. Używam "cameraman.tif".', target_file);
        filename = 'cameraman.tif';
        img_orig = imread(filename);
    end
    
    % Konwersja na GrayScale (wymuszenie 2D)
    if size(img_orig, 3) == 3
        img_orig = rgb2gray(img_orig);
        fprintf('Info: Konwersja RGB -> Gray dla "%s".\n', filename);
    end
    
    [rows, cols] = size(img_orig);
    num_pixels = rows * cols;
    
    % --- Generowanie 3 wersji modyfikacji (dla analizy wrażliwości) ---
    
    % 1. Start (Pierwszy piksel)
    img_mod_start = img_orig;
    img_mod_start(1) = bitxor(img_mod_start(1), 1);
    
    % 2. Mid (Środkowy piksel)
    img_mod_mid = img_orig;
    mid_idx = ceil(num_pixels / 2);
    img_mod_mid(mid_idx) = bitxor(img_mod_mid(mid_idx), 1);
    
    % 3. End (Ostatni piksel)
    img_mod_end = img_orig;
    img_mod_end(end) = bitxor(img_mod_end(end), 1);
    
    % --- OBLICZANIE STATYSTYK ORYGINAŁU (Tego brakowało) ---
    try
        % Zakładamy, że funkcje metryk są w ścieżce (dodane w main.m)
        entropy_orig = calculate_entropy(img_orig);
        [ch, cv, cd] = calculate_correlation(img_orig);
    catch ME
        warning('Nie można obliczyć metryk oryginału: %s', ME.message);
        entropy_orig = 0;
        ch = 0; cv = 0; cd = 0;
    end
    
    % --- Pakowanie danych do struktury ---
    data.img_orig = img_orig;
    data.size     = size(img_orig);
    data.filename = filename;
    data.img_flat = img_orig(:)';
    
    % Statystyki oryginału
    data.stats.entropy = entropy_orig;
    data.stats.corr    = [ch, cv, cd];
    
    % Wersje zmodyfikowane (Macierze i Płaskie)
    data.img_mod_start = img_mod_start;
    data.img_flat_mod_start = img_mod_start(:)';
    
    data.img_mod_mid = img_mod_mid;
    data.img_flat_mod_mid = img_mod_mid(:)';
    
    data.img_mod_end = img_mod_end;
    data.img_flat_mod_end = img_mod_end(:)';
    
    % Domyślne pola (dla kompatybilności wstecznej z resztą kodu)
    data.img_mod = img_mod_mid;       % Domyślnie bierzemy środkowy do wykresów
    data.img_flat_mod = img_mod_mid(:)';
end