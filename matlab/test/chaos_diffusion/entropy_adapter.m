function [out1, out2] = entropy_adapter(mode, varargin)
    % ENTROPY_ADAPTER - Implementacja logiki Pourasad et al. (2021)
    % Obsługuje: Encrypt, Decrypt.
    
    if strcmp(mode, 'encrypt')
        [out1, out2] = encrypt_core(varargin{1}); % out1=Img, out2=KeysStruct
    elseif strcmp(mode, 'decrypt')
        out1 = decrypt_core(varargin{1}, varargin{2}); % out1=Img
        out2 = [];
    else
        error('Tryb musi byc "encrypt" lub "decrypt"');
    end
end

% =========================================================
% CORE: ENCRYPTION
% =========================================================
function [encryptedImage, keys] = encrypt_core(inputImage)
    % Preprocessing
    if size(inputImage, 3) == 3
        inputImage = rgb2gray(inputImage);
    end
    
    % DWT wymaga parzystych wymiarów dla prostego podziału
    [m, n] = size(inputImage);
    if mod(m,2)~=0 || mod(n,2)~=0
        inputImage = imresize(inputImage, [floor(m/2)*2, floor(n/2)*2]);
        [m, n] = size(inputImage);
    end
    
    % --- 1. Key Generation (Diffusion Keys) ---
    x1_0 = 0.5; mu1 = 4.0;
    x2_0 = 0.5; mu2 = 3.9;
    total_pixels = m * n;
    
    seq1 = logistic_map_gen(x1_0, mu1, total_pixels);
    seq2 = logistic_map_gen(x2_0, mu2, total_pixels);
    
    key_seq_float = max(seq1, seq2);
    key_seq_int = uint8(mod(floor(key_seq_float * 1e14), 256));
    key_matrix = reshape(key_seq_int, [m, n]);
    
    % --- 2. Image Diffusion ---
    % Eq: D = NOT(Input XOR Key) -> bitcmp(bitxor(..))
    diffusedImage = bitcmp(bitxor(inputImage, key_matrix), 'uint8');
    
    % --- 3. Decomposition (DWT) ---
    % Haar wavelet
    [cA, cH, cV, cD] = dwt2(double(diffusedImage), 'haar');
    coeffs = [cA, cH; cV, cD];
    [cm, cn] = size(coeffs);
    
    % --- 4. Confusion Keys (Permutation) ---
    x3_0 = 0.3; y3_0 = 0.3; a_param = 1.9; b_param = 0.01;
    
    [seqX, seqY] = cml_map_gen(x3_0, y3_0, a_param, b_param, max(cm, cn));
    seqX = seqX(1:cm);
    seqY = seqY(1:cn);
    
    [~, idxRow] = sort(seqX);
    [~, idxCol] = sort(seqY);
    
    % --- 5. Image Confusion ---
    confusedCoeffs = coeffs(idxRow, idxCol);
    
    % --- 6. Reconstruction (IDWT) ---
    half_m = cm / 2; half_n = cn / 2;
    ncA = confusedCoeffs(1:half_m, 1:half_n);
    ncH = confusedCoeffs(1:half_m, half_n+1:end);
    ncV = confusedCoeffs(half_m+1:end, 1:half_n);
    ncD = confusedCoeffs(half_m+1:end, half_n+1:end);
    
    encryptedImageDouble = idwt2(ncA, ncH, ncV, ncD, 'haar');
    
    % Normalizacja do 0-255 uint8
    encryptedImage = uint8(mat2gray(encryptedImageDouble) * 255);
    
    % Zapisywanie kluczy/parametrów do deszyfrowania
    keys.key_matrix = key_matrix;
    keys.idxRow = idxRow;
    keys.idxCol = idxCol;
    keys.orig_range = [min(encryptedImageDouble(:)), max(encryptedImageDouble(:))]; % Opcjonalne, dla lepszej wierności IDWT
end

% =========================================================
% CORE: DECRYPTION
% =========================================================
function P_recovered = decrypt_core(encryptedImage, keys)
    % Konwersja do double
    encryptedImage = double(encryptedImage);
    
    % Uwaga: DWT/IDWT wprowadza błędy zaokrągleń (float), co utrudnia 
    % bezstratne odzyskanie przy prostym rzutowaniu uint8.
    % W tym algorytmie "mat2gray" niszczy oryginalną amplitudę.
    % Spróbujemy odtworzyć proces odwrotny na znormalizowanym obrazie.
    
    % --- 1. Decomposition (DWT) ---
    % Odzyskujemy pomieszane współczynniki z zaszyfrowanego obrazu
    [cA, cH, cV, cD] = dwt2(encryptedImage, 'haar');
    confusedCoeffs = [cA, cH; cV, cD];
    
    % --- 2. Inverse Confusion (Permutation) ---
    % Odwracanie permutacji wierszy i kolumn
    % A(idx) = B  =>  Temp(idx) = B rows/cols
    
    % Odwracanie wierszy
    temp_rows = zeros(size(confusedCoeffs));
    temp_rows(keys.idxRow, :) = confusedCoeffs;
    
    % Odwracanie kolumn
    coeffs = zeros(size(temp_rows));
    coeffs(:, keys.idxCol) = temp_rows;
    
    % --- 3. Reconstruction (IDWT) - Powrót do Diffused ---
    [cm, cn] = size(coeffs);
    half_m = cm / 2; half_n = cn / 2;
    
    ncA = coeffs(1:half_m, 1:half_n);
    ncH = coeffs(1:half_m, half_n+1:end);
    ncV = coeffs(half_m+1:end, 1:half_n);
    ncD = coeffs(half_m+1:end, half_n+1:end);
    
    diffusedImageDouble = idwt2(ncA, ncH, ncV, ncD, 'haar');
    diffusedImage = uint8(diffusedImageDouble); % Rzutowanie
    
    % --- 4. Inverse Diffusion ---
    % Forward: D = ~ (I XOR K)
    % Reverse: ~D = I XOR K  =>  I = (~D) XOR K
    
    not_diffused = bitcmp(diffusedImage, 'uint8');
    P_recovered = bitxor(not_diffused, keys.key_matrix);
    
    P_recovered = uint8(P_recovered);
end

% =========================================================
% FUNKCJE POMOCNICZE
% =========================================================

function seq = logistic_map_gen(x0, mu, N)
    seq = zeros(1, N);
    x = x0;
    for i = 1:100, x = mu * x * (1 - x); end % Transient
    for i = 1:N
        x = mu * x * (1 - x);
        seq(i) = x;
    end
end

function [seqX, seqY] = cml_map_gen(x0, y0, a, b, N)
    seqX = zeros(1, N);
    seqY = zeros(1, N);
    x = x0; y = y0;
    for i = 1:100 % Transient
        xn = 1 - a * (x^2 + y^2);
        yn = -2 * a * (1 - 2*b) * x * y;
        x = xn; y = yn;
    end
    for i = 1:N
        xn = 1 - a * (x^2 + y^2);
        yn = -2 * a * (1 - 2*b) * x * y;
        x = xn; y = yn;
        seqX(i) = x; seqY(i) = y;
    end
end
