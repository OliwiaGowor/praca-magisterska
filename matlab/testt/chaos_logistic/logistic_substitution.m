function img_out = logistic_substitution_encrypt(img_in, x1, x2, lam1, lam2)
    [M, N] = size(img_in);
    L = 256; % Grayscale levels
    
    % Generate Chaotic Matrices g1 and g2 [cite: 240]
    % Note: Discretizing chaotic sequence to be usable in modular arithmetic
    seq1 = two_stage_logistic(x1, lam1, M*N);
    seq2 = two_stage_logistic(x2, lam2, M*N);
    
    % Scale to integers (standard practice for chaotic image encryption)
    g1 = reshape(floor(seq1 * L), M, N);
    g2 = reshape(floor(seq2 * L), M, N);
    
    img_out = zeros(M, N);
    
    % Apply Formula (3) 
    % Note: Paper uses 0-based indexing for i,j formulas.
    for i = 1:M
        for j = 1:N
            % i_idx, j_idx are 0-based versions of i,j
            i_idx = i - 1;
            j_idx = j - 1;
            
            val = img_in(i,j) + (i_idx * g1(i,j)) + (j_idx * g2(i,j));
            img_out(i,j) = mod(val, L);
        end
    end
end

function img_out = logistic_substitution_decrypt(img_in, x1, x2, lam1, lam2)
    [M, N] = size(img_in);
    L = 256;
    
    % Re-generate identical Chaotic Matrices
    seq1 = two_stage_logistic(x1, lam1, M*N);
    seq2 = two_stage_logistic(x2, lam2, M*N);
    
    g1 = reshape(floor(seq1 * L), M, N);
    g2 = reshape(floor(seq2 * L), M, N);
    
    img_out = zeros(M, N);
    
    % Inverse operation of Formula (3)
    for i = 1:M
        for j = 1:N
            i_idx = i - 1;
            j_idx = j - 1;
            
            % Subtract the chaotic terms
            val = img_in(i,j) - (i_idx * g1(i,j)) - (j_idx * g2(i,j));
            
            % MATLAB mod() handles negatives correctly for this purpose 
            % (e.g., mod(-5, 256) -> 251)
            img_out(i,j) = mod(val, L);
        end
    end
end
