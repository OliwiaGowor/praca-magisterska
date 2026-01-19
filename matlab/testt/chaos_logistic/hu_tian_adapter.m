function [out1, out2] = hu_tian_adapter(mode, varargin)
    % HU_TIAN_ADAPTER - Implementation of Hu & Tian Algorithm (2020)
    % Integrates: logistic_substitution, two_stage_logistic, m_sequence_permute
    
    if strcmp(mode, 'encrypt')
        [out1, out2] = encrypt_core(varargin{1}); % out1=Img, out2=Keys
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
function [C, keys] = encrypt_core(P)
    % Parameters as defined in the provided main configuration
    keys.x1 = 0.41;
    keys.x2 = 0.87;
    keys.lambda1 = 3.95;
    keys.lambda2 = 3.80;
    keys.m_seed = 10;
    keys.iterations = 1;

    % Preprocessing (Grayscale)
    if size(P, 3) == 3
        P = rgb2gray(P);
    end
    P = double(P);
    
    C = P;
    
    % Main Encryption Loop
    for k = 1:keys.iterations
        % Step A: Substitution
        C = logistic_substitution(C, keys, 'encrypt');
        
        % Step B: Permutation (M-Sequence/Scrambling)
        C = m_sequence_permute(C, keys.m_seed, 'encrypt');
    end
    
    C = uint8(C);
end

% =========================================================
% CORE: DECRYPTION
% =========================================================
function P_recovered = decrypt_core(C, keys)
    C = double(C);
    P_recovered = C;
    
    % Reverse Decryption Loop
    for k = 1:keys.iterations
        % Reverse Step B: Permutation
        P_recovered = m_sequence_permute(P_recovered, keys.m_seed, 'decrypt');
        
        % Reverse Step A: Substitution
        P_recovered = logistic_substitution(P_recovered, keys, 'decrypt');
    end
    
    P_recovered = uint8(P_recovered);
end

% =========================================================
% HELPER FUNCTIONS
% =========================================================

function seq = two_stage_logistic(x0, gamma, len)
    % Generates chaotic sequence using Two-Stage Logistic Map
    seq = zeros(len, 1);
    x = x0;
    for k = 1:len
        if x >= 0 && x < 0.5
            x = 4 * gamma * x * (0.5 - x);
        elseif x >= 0.5 && x <= 1
            x = 1 - 4 * gamma * x * (x - 0.5);
        end
        seq(k) = x;
    end
end

function img_out = logistic_substitution(img_in, keys, mode)
    % Integrated substitution function
    [M, N] = size(img_in);
    L = 256;
    
    % Generate chaotic matrices
    seq1 = two_stage_logistic(keys.x1, keys.lambda1, M*N);
    seq2 = two_stage_logistic(keys.x2, keys.lambda2, M*N);
    
    g1 = reshape(floor(seq1 * L), M, N);
    g2 = reshape(floor(seq2 * L), M, N);
    
    img_out = zeros(M, N);
    
    if strcmp(mode, 'encrypt')
        for i = 1:M
            for j = 1:N
                i_idx = i - 1; j_idx = j - 1;
                val = img_in(i,j) + (i_idx * g1(i,j)) + (j_idx * g2(i,j));
                img_out(i,j) = mod(val, L);
            end
        end
    else % decrypt
        for i = 1:M
            for j = 1:N
                i_idx = i - 1; j_idx = j - 1;
                val = img_in(i,j) - (i_idx * g1(i,j)) - (j_idx * g2(i,j));
                img_out(i,j) = mod(val, L);
            end
        end
    end
end

function img_out = m_sequence_permute(img_in, seed_val, mode)
    % Position scrambling based on M-Sequence (simulated via randperm)
    [M, N] = size(img_in);
    num_pixels = M * N;
    
    % Deterministic permutation based on seed
    s = rng; % Save current generator state
    rng(seed_val); 
    perm_order = randperm(num_pixels);
    rng(s); % Restore state (to preserve randomness elsewhere)
    
    img_linear = img_in(:);
    img_out_linear = zeros(num_pixels, 1);
    
    if strcmp(mode, 'encrypt')
        img_out_linear(perm_order) = img_linear;
    elseif strcmp(mode, 'decrypt')
        img_out_linear = img_linear(perm_order);
    end
    
    img_out = reshape(img_out_linear, M, N);
end
