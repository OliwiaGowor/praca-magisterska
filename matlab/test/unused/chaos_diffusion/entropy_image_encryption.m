function entropy_chaos_implementation()
    % Implementation of "A New Algorithm for Digital Image Encryption Based on Chaos Theory"
    % Source: Entropy 2021, 23, 341. [cite: 7]
    
    clc; clear; close all;

    %% 1. Input Image Loading [cite: 195]
    % "In the first step, a grayscale image G is arranged."
    try
        G = imread('cameraman.tif'); % Load standard test image
    catch
        error('Image not found. Please add "cameraman.tif" or change filename.');
    end
    
    if size(G, 3) == 3
        G = rgb2gray(G);
    end
    G = imresize(G, [256, 256]); % Resize to m x n
    [m, n] = size(G);
    
    figure('Name', 'Entropy Algorithm Steps', 'Color', 'w');
    subplot(2,2,1); imshow(G); title('Original Image (G)');
    
    %% Step NO.1: Image Diffusion (Logistic Maps + XNOR) 
    % "By evaluating two logistic maps, a chaotic sequence is generated."
    
    % Parameters from Table 2 [cite: 284-290]
    x1_0 = 0.5;   mu1 = 4.0;
    x2_0 = 0.5;   mu2 = 3.9;
    
    % Generate chaotic sequences (Eq 1 & 2)
    seq_len = m * n;
    X1 = logistic_map(x1_0, mu1, seq_len);
    X2 = logistic_map(x2_0, mu2, seq_len);
    
    % Combine sequences to create Secure Key
    % The paper implies comparing and selecting values (Section 2.2 Step 2)
    Key_Seq = max(X1, X2); 
    
    % Quantize to uint8 (0-255) for XNOR operation
    Key_Int = uint8(mod(floor(Key_Seq * 1e14), 256));
    Key_Matrix = reshape(Key_Int, [m, n]);
    
    % "Making XNOR with the primary image, the diffusion is terminated." [cite: 196-197]
    % XNOR is ~(A XOR B). In MATLAB: bitcmp(bitxor(A, B))
    Diffused_Image = bitcmp(bitxor(G, Key_Matrix));

    %% Step NO.2: Wavelet Decomposition 
    % "wavelet decomposition is performed... coefficient is extracted, registered as ca1"
    
    % Using Haar wavelet (standard for simple DWT)
    [cA, cH, cV, cD] = dwt2(double(Diffused_Image), 'haar');
    
    % Combine coefficients into one matrix 'R' for shuffling
    % This creates the matrix to be confused (registered as 'ca1' in text)
    ca1 = [cA, cH; cV, cD]; 
    [mc, nc] = size(ca1);
    
    %% Step NO.3: Position Confusion (CML Map) 
    % "Utilizing a two-dimensional hyper-chaotic map CML... position confusion is performed"
    
    % Initial values from Table 2 [cite: 284-290]
    x3_0 = 0.3;
    y3_0 = 0.3;
    % Parameters a and b (Standard CML values as they are not in Table 2)
    a = 1.9; b = 0.01; 
    
    % Generate chaotic coordinate sequences (Eq 3 & 4)
    [Seq_X, Seq_Y] = cml_map(x3_0, y3_0, a, b, max(mc, nc));
    
    % Use only required length
    Seq_Row = Seq_X(1:mc);
    Seq_Col = Seq_Y(1:nc);
    
    % Sort to get permutation indices (Step NO.2 in Section 2.3)
    [~, Row_Idx] = sort(Seq_Row);
    [~, Col_Idx] = sort(Seq_Col);
    
    % Shuffle the coefficients 'ca1' (Eq 5)
    Confused_ca1 = ca1(Row_Idx, Col_Idx);

    subplot(2,2,2); imshow(Confused_ca1, []); title('Confused Coefficients (ca1)');

    %% Step NO.4: Reconstruction 
    % "confused image can be rebuilt by wavelet. After all, the encrypted image is obtained."
    
    % Separate the shuffled coefficients back into sub-bands
    half_m = mc / 2;
    half_n = nc / 2;
    
    new_cA = Confused_ca1(1:half_m, 1:half_n);
    new_cH = Confused_ca1(1:half_m, half_n+1:end);
    new_cV = Confused_ca1(half_m+1:end, 1:half_n);
    new_cD = Confused_ca1(half_m+1:end, half_n+1:end);
    
    % Inverse Discrete Wavelet Transform
    Encrypted_Double = idwt2(new_cA, new_cH, new_cV, new_cD, 'haar');
    
    % Normalize and cast to uint8
    Encrypted_Image = uint8(mat2gray(Encrypted_Double) * 255);
    
    subplot(2,2,3); imshow(Encrypted_Image); title('Final Encrypted Image');
    subplot(2,2,4); imhist(Encrypted_Image); title('Histogram');

    %% Analysis Results
    fprintf('Entropy of Original: %.4f\n', entropy(G));
    fprintf('Entropy of Encrypted: %.4f\n', entropy(Encrypted_Image));
    fprintf('NPCR (Original vs Encrypted): %.4f%%\n', calc_npcr(G, Encrypted_Image));
end

%% Helper Functions

function seq = logistic_map(x0, mu, N)
    % Logistic Map Equation (1) & (2) [cite: 118-120]
    seq = zeros(1, N);
    x = x0;
    % Discard transients
    for i=1:100
        x = mu * x * (1 - x);
    end
    % Generate
    for i=1:N
        x = mu * x * (1 - x);
        seq(i) = x;
    end
end

function [seqX, seqY] = cml_map(x0, y0, a, b, N)
    % CML Map Equations (3) & (4) [cite: 149-150]
    seqX = zeros(1, N);
    seqY = zeros(1, N);
    x = x0; y = y0;
    % Discard transients
    for i=1:100
        xn = 1 - a*(x^2 + y^2);
        yn = -2*a*(1 - 2*b)*x*y;
        x = xn; y = yn;
    end
    % Generate
    for i=1:N
        xn = 1 - a*(x^2 + y^2);
        yn = -2*a*(1 - 2*b)*x*y;
        x = xn; y = yn;
        seqX(i) = x;
        seqY(i) = y;
    end
end

function val = calc_npcr(img1, img2)
    % NPCR Calculation [cite: 240-242]
    d = img1 ~= img2;
    val = sum(d(:)) / numel(img1) * 100;
end
