function [psnr_val, ssim_val, sd_val] = calculate_alteration_metrics(img_original, img_decrypted)
%
% Oblicza metryki jakości deszyfrowania po ataku na szyfrogram.
% Porównuje oryginalny obraz A z odszyfrowanym obrazem D.
% WYMAGA: Image Processing Toolbox (dla psnr i ssim)
%
% img_original (A): Oryginalny obraz jawny.
% img_decrypted (D): Obraz odszyfrowany z uszkodzonego szyfrogramu.
%

    % Upewnij się, że oba obrazy są uint8 do porównania
    img_original = uint8(img_original);
    img_decrypted = uint8(img_decrypted);

    % 1. PSNR (Peak Signal-to-Noise Ratio) - Im wyżej, tym lepiej.
    % Formuła (25) jest zaimplementowana w funkcji psnr()
    psnr_val = psnr(img_decrypted, img_original);

    % 2. SSIM (Structural Similarity Index) - Im bliżej 1, tym lepiej.
    % Formuła (24) jest zaimplementowana w funkcji ssim()
    ssim_val = ssim(img_decrypted, img_original);

    % 3. SD (Spectral Distortion) - Im niżej, tym lepiej.
    % Formuła (26)
    F_A = fft2(double(img_original));
    F_D = fft2(double(img_decrypted));
    
    % 1/(M*N) * suma(...) to po prostu średnia 'mean'
    sd_val = mean(abs(F_A - F_D), 'all');

end