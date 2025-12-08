function [t_enc, t_dec, C1, C2, PT] = algo_hu_tian(data)
    % ALGO_HU_TIAN - Wrapper dla algorytmu Hu & Tian (Two-Stage Logistic + M-Sequence)
    % Zgodny ze standardem testowym.
    
    % 1. Przygotowanie danych
    img_orig = data.img_orig;
    img_mod  = data.img_mod;
    
    % 2. Szyfrowanie (Oryginał) - pomiar czasu
    tic;
    % Adapter zwraca zaszyfrowany obraz (ct1) i strukturę kluczy (keys)
    [ct1, keys] = hu_tian_adapter('encrypt', img_orig);
    t_enc = toc;
    
    % 3. Szyfrowanie (Zmodyfikowany) - dla NPCR/UACI
    % Używamy tych samych parametrów/kluczy generowanych wewnątrz adaptera
    % (ale dla poprawności testu NPCR, algorytm powinien generować strumień deterministycznie)
    [ct2, ~] = hu_tian_adapter('encrypt', img_mod);
    
    % 4. Deszyfrowanie - pomiar czasu
    tic;
    pt_img = hu_tian_adapter('decrypt', ct1, keys);
    t_dec = toc;
    
    % 5. Formatowanie wyników
    C1 = uint8(ct1);
    C2 = uint8(ct2);
    PT = uint8(pt_img);
end
