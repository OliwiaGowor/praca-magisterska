function [t_enc, t_dec, C1, C2, PT] = algo_entropy_pourasad(data)
    % ALGO_ENTROPY_POURASAD - Wrapper dla algorytmu Pourasad et al. (2021)
    % Bazuje na DWT, mapach chaotycznych (Logistic + CML) i operacjach bitowych.
    
    % 1. Przygotowanie danych
    img_orig = data.img_orig;
    img_mod  = data.img_mod;
    
    % Algorytm pracuje na Grayscale (wymóg DWT i logiki autora).
    % Jeśli wejście jest RGB, konwertujemy (wynik PT będzie grayscale).
    
    % 2. Szyfrowanie (Oryginał)
    tic;
    [ct1, keys] = entropy_adapter('encrypt', img_orig);
    t_enc = toc;
    
    % 3. Szyfrowanie (Zmodyfikowany) - dla metryk NPCR/UACI
    % Używamy tych samych kluczy/parametrów, aby test był miarodajny
    [ct2, ~] = entropy_adapter('encrypt', img_mod);
    
    % 4. Deszyfrowanie
    tic;
    pt_img = entropy_adapter('decrypt', ct1, keys);
    t_dec = toc;
    
    % 5. Formatowanie wyników
    C1 = uint8(ct1);
    C2 = uint8(ct2);
    PT = uint8(pt_img);
end
