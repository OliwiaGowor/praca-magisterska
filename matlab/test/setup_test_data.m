function data = setup_test_data(custom_img)
    % SETUP_TEST_DATA - Przygotowuje obrazy do testów
    
    if nargin < 1 || isempty(custom_img)
        try
            img_orig = imread('cameraman.tif');
        catch
            img_orig = uint8(randi([0 255], 256, 256));
        end
    else
        img_orig = custom_img;
    end

    if size(img_orig, 3) == 3
        img_orig = rgb2gray(img_orig);
    end
    
    img_orig = uint8(img_orig);
    
    % Obraz modyfikowany (1 bit)
    img_mod = img_orig;
    img_mod(end) = bitxor(img_mod(end), 1);
    
    % --- WAŻNE: Wersje 1D dla funkcji MEX ---
    img_flat = img_orig(:)';
    img_flat_mod = img_mod(:)';
    
    % Analysis of original
    entropy_orig = calculate_entropy(img_orig);
    [ch, cv, cd] = calculate_correlation(img_orig);
    
    % Pack into output structure
    data.img_orig = img_orig;
    data.img_mod = img_mod;
    data.stats.entropy = entropy_orig;
    data.stats.corr = [ch, cv, cd];
    data.img_flat = img_flat;          % <--- To pole jest wymagane
    data.img_flat_mod = img_flat_mod;  % <--- To pole jest wymagane
    data.size = size(img_orig);
end