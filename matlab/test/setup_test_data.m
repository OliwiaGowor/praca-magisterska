function data = setup_test_data()
    % Loads image, creates modified version, and calculates original statistics
    
    try
        img_orig = imread('cameraman.tif');
        % img_orig = imresize(img_orig, [64 64]); % Optional resizing
    catch
        warning('Brak pliku cameraman.tif, generowanie szumu.');
        img_orig = uint8(randi([0 255], 256, 256));
    end
    
    img_size = size(img_orig);
    
    % --- Modified image (for NPCR/UACI) ---
    img_mod = img_orig;
    img_mod(1) = bitxor(img_mod(1), 1); % Change 1 bit
    
    % Flattened versions (1D) for C/MEX functions
    img_flat = img_orig(:)';
    img_flat_mod = img_mod(:)';
    
    % Analysis of original
    entropy_orig = calculate_entropy(img_orig);
    [ch, cv, cd] = calculate_correlation(img_orig);
    
    % Pack into output structure
    data.img_orig = img_orig;
    data.img_mod = img_mod;
    data.img_flat = img_flat;
    data.img_flat_mod = img_flat_mod;
    data.size = img_size;
    data.stats.entropy = entropy_orig;
    data.stats.corr = [ch, cv, cd];
end
