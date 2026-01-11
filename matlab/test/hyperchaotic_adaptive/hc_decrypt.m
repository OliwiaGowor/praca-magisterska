function plain_img = hc_decrypt(cipher, chaotic_seq, sz)
    % HC_DECRYPT Główna procedura deszyfrowania: Odwrotna Konfuzja -> Odwrotna Dyfuzja
    
    rows = sz(1); cols = sz(2); 
    if length(sz) > 2, chans = sz(3); else, chans = 1; end
    num_pixels = rows * cols;
    
    plain_img = zeros(sz, 'uint8');
    seq_idx = 1;
    
    for k = 1:chans
        channel = double(cipher(:,:,k));
        
        % Pobieranie tych samych fragmentów sekwencji co przy szyfrowaniu
        seq_hist = chaotic_seq(seq_idx : seq_idx + 255);
        seq_idx = seq_idx + 256;
        
        seq_conf = chaotic_seq(seq_idx : seq_idx + num_pixels - 1);
        seq_idx = seq_idx + num_pixels;
        
        % --- KROK 1: ODWROTNA KONFUZJA ---
        % Wejście do descramble musi być spłaszczone
        descrambled_flat = hc_confusion_inv(channel(:), seq_conf);
        
        % --- KROK 2: ODWROTNA DYFUZJA ---
        restored_flat = hc_diffusion_inv(descrambled_flat, seq_hist);
        
        plain_img(:,:,k) = uint8(reshape(restored_flat, rows, cols));
    end
end