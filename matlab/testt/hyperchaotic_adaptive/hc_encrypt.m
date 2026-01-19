function cipher_img = hc_encrypt(img, chaotic_seq)
    % HC_ENCRYPT Główna procedura szyfrowania: Dyfuzja -> Konfuzja
    
    [rows, cols, chans] = size(img);
    num_pixels = rows * cols;
    cipher_img = zeros(size(img), 'uint8');
    
    seq_idx = 1;
    
    for k = 1:chans
        channel = double(img(:,:,k));
        
        % --- KROK 1: SELF-ADAPTIVE DIFFUSION ---
        % Pobierz 256 wartości do inicjalizacji histogramu
        seq_hist = chaotic_seq(seq_idx : seq_idx + 255);
        seq_idx = seq_idx + 256;
        
        diffused_flat = hc_diffusion(channel, seq_hist);
        
        % --- KROK 2: CONFUSION ---
        % Pobierz num_pixels wartości do permutacji
        seq_conf = chaotic_seq(seq_idx : seq_idx + num_pixels - 1);
        seq_idx = seq_idx + num_pixels;
        
        confused_flat = hc_confusion(diffused_flat, seq_conf);
        
        cipher_img(:,:,k) = uint8(reshape(confused_flat, rows, cols));
    end
end