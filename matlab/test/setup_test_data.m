function data = setup_test_data()
    % SETUP_TEST_DATA - Prepares images (Original + 3 Modifications) + Statistics
    
    target_file = 'images/4.2.03.tiff'; 
    
    if exist(target_file, 'file')
        filename = target_file;
        img_orig = imread(filename);
    else
        warning('Plik %s nie istnieje. Używam "cameraman.tif".', target_file);
        filename = 'cameraman.tif';
        img_orig = imread(filename);
    end
    
    % Conversion to GrayScale (forcing 2D)
    if size(img_orig, 3) == 3
        img_orig = rgb2gray(img_orig);
        fprintf('Info: Konwersja RGB -> Gray dla "%s".\n', filename);
    end
    
    [rows, cols] = size(img_orig);
    num_pixels = rows * cols;
    
    % --- Generating 3 modification versions (for sensitivity analysis) ---
    
    % 1. Start (First pixel)
    img_mod_start = img_orig;
    img_mod_start(1) = bitxor(img_mod_start(1), 1);
    
    % 2. Mid (Middle pixel)
    img_mod_mid = img_orig;
    mid_idx = ceil(num_pixels / 2);
    img_mod_mid(mid_idx) = bitxor(img_mod_mid(mid_idx), 1);
    
    % 3. End (Last pixel)
    img_mod_end = img_orig;
    img_mod_end(end) = bitxor(img_mod_end(end), 1);
    
    % --- CALCULATING ORIGINAL STATISTICS ---
    try
        entropy_orig = calculate_entropy(img_orig);
        [ch, cv, cd] = calculate_correlation(img_orig);
    catch ME
        warning('Nie można obliczyć metryk oryginału: %s', ME.message);
        entropy_orig = 0;
        ch = 0; cv = 0; cd = 0;
    end
    
    % --- Packing data into structure ---
    data.img_orig = img_orig;
    data.size     = size(img_orig);
    data.filename = filename;
    data.img_flat = img_orig(:)';
    
    % Original statistics
    data.stats.entropy = entropy_orig;
    data.stats.corr    = [ch, cv, cd];
    
    % Modified versions (Matrices and Flat)
    data.img_mod_start = img_mod_start;
    data.img_flat_mod_start = img_mod_start(:)';
    
    data.img_mod_mid = img_mod_mid;
    data.img_flat_mod_mid = img_mod_mid(:)';
    
    data.img_mod_end = img_mod_end;
    data.img_flat_mod_end = img_mod_end(:)';
    
    data.img_mod = img_mod_mid;
    data.img_flat_mod = img_mod_mid(:)';
end
