function [psnr_val, ssim_val, sd_val] = calculate_alteration_metrics(img_original, img_decrypted)
%
% Calculates decryption quality metrics after a ciphertext attack.
% Compares original image A with decrypted image D.
% REQUIRES: Image Processing Toolbox (for psnr and ssim)
%
% img_original (A): Original plaintext image.
% img_decrypted (D): Decrypted image from corrupted ciphertext.
%

    % Ensure both images are uint8 for comparison
    img_original = uint8(img_original);
    img_decrypted = uint8(img_decrypted);

    % 1. PSNR (Peak Signal-to-Noise Ratio) - Higher is better.
    psnr_val = psnr(img_decrypted, img_original);

    % 2. SSIM (Structural Similarity Index) - Closer to 1 is better.
    ssim_val = ssim(img_decrypted, img_original);

    % 3. SD (Spectral Distortion) - Lower is better.
    F_A = fft2(double(img_original));
    F_D = fft2(double(img_decrypted));
    
    % Calculate the mean absolute difference in the frequency domain
    sd_val = mean(abs(F_A - F_D), 'all');

end
