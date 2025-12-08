function img_out = m_sequence_permute(img_in, seed_val, mode)
    % Implements position scrambling (permutation).
    % Section 3.2.1  describes using M-sequence ergodicity 
    % to shuffle positions.
    
    [M, N] = size(img_in);
    num_pixels = M * N;
    
    % Generate permutation vector
    % We use a seeded RNG to simulate the M-sequence generator's deterministic nature
    rng(seed_val); 
    perm_order = randperm(num_pixels);
    
    img_linear = img_in(:);
    img_out_linear = zeros(num_pixels, 1);
    
    if strcmp(mode, 'encrypt')
        % Move pixels TO the permuted positions
        % img_out(p(i)) = img_in(i)
        % (Forward scrambling)
        img_out_linear(perm_order) = img_linear;
        
    elseif strcmp(mode, 'decrypt')
        % Move pixels BACK from the permuted positions
        % img_out(i) = img_in(p(i))
        img_out_linear = img_linear(perm_order);
    end
    
    img_out = reshape(img_out_linear, M, N);
end
