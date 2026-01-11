function [t_enc, t_dec, C1, C2, PT] = algo_des_c(data)
    % ALGO_DES_C - Wrapper dla DES-CBC (OpenSSL via MEX)
    % Wymaga skompilowanego crypto_mex z obsługą 'legacy' (OpenSSL 3.0+)
    
    % 1. Przygotowanie danych (spłaszczenie do 1D)
    if isfield(data, 'img_flat')
        p_flat = data.img_flat;
        p_mod_flat = data.img_flat_mod;
    else
        p_flat = data.img_orig(:)';
        p_mod_flat = data.img_mod(:)';
    end
    
    sz = data.size;
    
    % DES: Klucz 8 bajtów (64 bity, z czego 56 efektywnych), IV 8 bajtów
    key = uint8(randi([0 255], 1, 8));
    iv  = uint8(randi([0 255], 1, 8));
    
    algo_name = 'des-cbc';
    
    % 2. Szyfrowanie (Oryginał)
    tic;
    ct1_flat = crypto_mex(p_flat, key, iv, algo_name, true);
    t_enc = toc;
    
    % 3. Szyfrowanie (Zmodyfikowany)
    ct2_flat = crypto_mex(p_mod_flat, key, iv, algo_name, true);
    
    % 4. Deszyfrowanie
    tic;
    pt_flat = crypto_mex(ct1_flat, key, iv, algo_name, false);
    t_dec = toc;
    
    % 5. Formatowanie wyników
    % DES (blokowy) dodaje padding (PKCS7). Wynik jest dłuższy niż oryginał.
    % Do wizualizacji i metryk przycinamy go do rozmiaru obrazu.
    
    len = prod(sz);
    
    % Zabezpieczenie przed błędami wymiarów
    if numel(ct1_flat) >= len
        C1 = reshape(ct1_flat(1:len), sz);
    else
        % Fallback w razie błędu (np. pusty wynik)
        C1 = zeros(sz, 'uint8');
    end
    
    if numel(ct2_flat) >= len
        C2 = reshape(ct2_flat(1:len), sz);
    else
        C2 = zeros(sz, 'uint8');
    end
    
    % Odszyfrowany obraz powinien być idealnie równy oryginałowi (po zdjęciu paddingu)
    if numel(pt_flat) >= len
        PT = reshape(pt_flat(1:len), sz);
    else
        PT = zeros(sz, 'uint8');
    end
end