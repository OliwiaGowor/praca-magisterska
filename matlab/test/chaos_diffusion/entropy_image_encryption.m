function entropy_image_encryption()
% ENTROPY_IMAGE_ENCRYPTION Implements the algorithm from Pourasad et al. (2021)
% "A New Algorithm for Digital Image Encryption Based on Chaos Theory"
% Source: Entropy 2021, 23, 341.
%
% Steps:
% 1. Image Diffusion (Logistic Map + XNOR)
% 2. Discrete Wavelet Transform (DWT)
% 3. Confusion (2D CML Map Permutation on Coefficients)
% 4. Inverse Discrete Wavelet Transform (IDWT)

    clc; clear; close all;

    %% 1. Input Image Loading
    % Load a standard image (e.g., Cameraman or Lena)
    try
        inputImage = imread('cameraman.tif'); % Built-in MATLAB image
    catch
        error('Image not found. Please ensure "cameraman.tif" or another image is in path.');
    end
    
    % Resize to even dimensions for simpler DWT handling
    [m, n, ch] = size(inputImage);
    if ch > 1
        inputImage = rgb2gray(inputImage);
    end
    inputImage = imresize(inputImage, [256, 256]);
    [m, n] = size(inputImage);
    
    fprintf('Processing Image Size: %dx%d\n', m, n);

    %% 2. Image Diffusion (Section 2.2)
    % Parameters from Table 2
    x1_0 = 0.5;   mu1 = 4.0;
    x2_0 = 0.5;   mu2 = 3.9;
    
    % Generate two chaotic sequences using Logistic Map
    total_pixels = m * n;
    seq1 = logistic_map_gen(x1_0, mu1, total_pixels);
    seq2 = logistic_map_gen(x2_0, mu2, total_pixels);
    
    % Compare and choose the larger value (Step 2)
    key_seq_float = max(seq1, seq2);
    
    % Convert float sequence to uint8 integers for XNOR operation
    % Assumption: floor(val * 1e14) mod 256 (Common practice in chaos encryption)
    key_seq_int = uint8(mod(floor(key_seq_float * 1e14), 256));
    key_matrix = reshape(key_seq_int, [m, n]);
    
    % Perform XNOR Operation (Step 3)
    % XNOR is ~(A XOR B). In MATLAB for uint8: bitcmp(bitxor(A, B))
    diffusedImage = bitcmp(bitxor(inputImage, key_matrix), 'uint8');
    
    %% 3. Image Decomposition (DWT) (Section 3 Step 2)
    % Use Haar wavelet for decomposition
    [cA, cH, cV, cD] = dwt2(double(diffusedImage), 'haar');
    
    % Construct the coefficient matrix (Concatenate sub-bands)
    % Figure 2 structure: [LL, LH; HL, HH] -> [cA, cV; cH, cD] approx
    % Standard MATLAB arrangement: [cA, cH; cV, cD]
    coeffs = [cA, cH; cV, cD];
    [cm, cn] = size(coeffs);
    
    %% 4. Image Confusion (Section 2.3 & 3 Step 3)
    % Using 2D Hyper-Chaotic Map (CML)
    % Initial values from Table 2
    x3_0 = 0.3;
    y3_0 = 0.3;
    
    % Constants a and b (Not in table, assumed standard chaotic values)
    a_param = 1.9; 
    b_param = 0.01;
    
    % Generate chaotic sequences X and Y
    % Lengths must match the coefficient matrix dimensions
    [seqX, seqY] = cml_map_gen(x3_0, y3_0, a_param, b_param, max(cm, cn));
    
    seqX = seqX(1:cm); % Row indices
    seqY = seqY(1:cn); % Col indices
    
    % Sort sequences to get permutation indices
    [~, idxRow] = sort(seqX);
    [~, idxCol] = sort(seqY);
    
    % Permute the coefficient matrix
    % R(i,j) = R(w2(i), w3(j))
    confusedCoeffs = coeffs(idxRow, idxCol);
    
    %% 5. Image Reconstruction (IDWT) (Section 3 Step 4)
    % Split confused coefficients back into sub-bands
    half_m = cm / 2;
    half_n = cn / 2;
    
    ncA = confusedCoeffs(1:half_m, 1:half_n);
    ncH = confusedCoeffs(1:half_m, half_n+1:end);
    ncV = confusedCoeffs(half_m+1:end, 1:half_n);
    ncD = confusedCoeffs(half_m+1:end, half_n+1:end);
    
    % Inverse DWT
    encryptedImageDouble = idwt2(ncA, ncH, ncV, ncD, 'haar');
    
    % Normalize/Cast to uint8 for display (Encrypted image looks like noise)
    encryptedImage = uint8(mat2gray(encryptedImageDouble) * 255);

    %% 6. Results and Histogram Analysis
    figure('Name', 'Entropy 2021 Algorithm Implementation', 'Color', 'w');
    
    subplot(2,3,1); imshow(inputImage); title('Original Image');
    subplot(2,3,2); imshow(diffusedImage); title('After Diffusion (Logistic)');
    subplot(2,3,3); imshow(encryptedImage); title('Final Encrypted (DWT+CML)');
    
    subplot(2,3,4); imhist(inputImage); title('Histogram: Original');
    subplot(2,3,5); imhist(diffusedImage); title('Histogram: Diffused');
    subplot(2,3,6); imhist(encryptedImage); title('Histogram: Encrypted');
    
    % Calculation of Entropy
    fprintf('Entropy of Original: %.4f\n', entropy(inputImage));
    fprintf('Entropy of Encrypted: %.4f\n', entropy(encryptedImage));
    
    % Calculation of NPCR (Number of Pixels Change Rate)
    % Comparing Diffused vs Encrypted (as an example of change)
    npcr_val = npcr_calc(inputImage, encryptedImage);
    fprintf('NPCR (Original vs Encrypted): %.4f%%\n', npcr_val);

end

%% Helper Functions

function seq = logistic_map_gen(x0, mu, N)
    % Eq (1) & (2): x(n+1) = mu * x(n) * (1 - x(n))
    seq = zeros(1, N);
    x = x0;
    % Discard transients (optional, but standard practice)
    for i = 1:100
        x = mu * x * (1 - x);
    end
    for i = 1:N
        x = mu * x * (1 - x);
        seq(i) = x;
    end
end

function [seqX, seqY] = cml_map_gen(x0, y0, a, b, N)
    % Eq (3) & (4): CML Map / TD-ERCS
    % x(n+1) = 1 - a*(x(n)^2 + y(n)^2)
    % y(n+1) = -2*a*(1 - 2*b)*x(n)*y(n)
    
    seqX = zeros(1, N);
    seqY = zeros(1, N);
    x = x0;
    y = y0;
    
    % Discard transients
    for i = 1:100
        xn = 1 - a * (x^2 + y^2);
        yn = -2 * a * (1 - 2*b) * x * y;
        x = xn;
        y = yn;
    end
    
    for i = 1:N
        xn = 1 - a * (x^2 + y^2);
        yn = -2 * a * (1 - 2*b) * x * y;
        x = xn;
        y = yn;
        
        seqX(i) = x;
        seqY(i) = y;
    end
end

function val = npcr_calc(img1, img2)
    % NPCR Calculation
    diff = img1 ~= img2;
    val = sum(diff(:)) / numel(img1) * 100;
end
