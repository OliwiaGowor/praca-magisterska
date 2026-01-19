function [t_enc, t_dec, C1, C2, PT] = algo_chaos_bit_shuffling_original(data)
    % ALGO_CHAOS_BIT_SHUFFLING_ORIGINAL - Wrapper dla wersji oryginalnej (wolnej)
    % Porównujemy ją z wersją _optimized.
    
    img_orig = data.img_orig;
    img_mod  = data.img_mod;
    
    % 1. Szyfrowanie (Oryginał)
    tic;
    % Wywołujemy plik 'chaos_bit_shuffling.m' (nie _optimized)
    [ct1, keys] = chaos_bit_shuffling('encrypt', img_orig);
    t_enc = toc;
    
    % 2. Szyfrowanie (Zmodyfikowany) - dla metryk
    [ct2, ~] = chaos_bit_shuffling('encrypt', img_mod);
    
    % 3. Deszyfrowanie
    tic;
    pt_img = chaos_bit_shuffling('decrypt', ct1, keys);
    t_dec = toc;
    
    % 4. Konwersja wyników
    % Oryginalny algorytm często zwraca double, rzutujemy na uint8
    C1 = uint8(ct1);
    C2 = uint8(ct2);
    PT = uint8(pt_img);
end
