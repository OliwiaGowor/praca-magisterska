function [t_enc, t_dec, C1, C2, PT] = algo_aes_pixel(data)
    % AES w trybie "Pikselowym" (ECB) - Wersja Zoptymalizowana (Wektorowa)
    
    % Generuj klucz 256-bit
    key = uint8(randi([0 255], 1, 32)); 
    
    input_flat = data.img_flat;     % Wektor 1D
    input_mod  = data.img_flat_mod; % Wektor 1D
    sz = data.size;
    
    % Inicjalizacja tablic (robimy to raz, poza pomiarem czasu szyfrowania per se,
    % lub wliczamy - zależy od metodologii. Tutaj wliczamy dla uczciwości wobec DES).
    [w, Nr, s_box, inv_s_box, m9, m11, m13, m14] = aes_init(key);
    
    % --- SZYFROWANIE (Wektoryzowane) ---
    tic;
    % AES wymaga danych będących wielokrotnością 16 bajtów. 
    % Paddingujemy cały wektor raz.
    [padded_input, pad_len] = pkcs7_pad(input_flat);
    
    % Wywołanie nowej, szybkiej funkcji
    ct1_padded = aes_cipher_vectorized(padded_input, w, Nr, s_box);
    t_enc = toc;
    
    % Szyfrowanie drugiego obrazu (bez pomiaru czasu)
    [padded_mod, ~] = pkcs7_pad(input_mod);
    ct2_padded = aes_cipher_vectorized(padded_mod, w, Nr, s_box);
    
    % --- DESZYFROWANIE (Wektoryzowane) ---
    tic;
    % Używamy istniejącej funkcji wektorowej (była używana w CBC, tu też zadziała idealnie)
    pt_padded = aes_inv_cipher_vectorized(ct1_padded, w, Nr, inv_s_box, m9, m11, m13, m14);
    
    % Usunięcie paddingu
    if pad_len > 0
        pt_flat = pt_padded(1:end-pad_len);
    else
        pt_flat = pt_padded;
    end
    t_dec = toc;
    
    % --- Formatowanie wyników ---
    % Przycinamy szyfrogram do rozmiaru obrazu dla celów wizualizacji (zgodnie z DES)
    num_pixels = prod(sz);
    C1 = reshape(ct1_padded(1:num_pixels), sz);
    C2 = reshape(ct2_padded(1:num_pixels), sz);
    
    % Odszyfrowany obraz
    PT = reshape(pt_flat(1:num_pixels), sz);
end

function [padded, p_len] = pkcs7_pad(data)
    % Pomocnicza funkcja do paddingu
    len = length(data);
    p_len = 16 - mod(len, 16);
    if p_len == 0, p_len = 16; end
    padded = [data, repmat(uint8(p_len), 1, p_len)];
end