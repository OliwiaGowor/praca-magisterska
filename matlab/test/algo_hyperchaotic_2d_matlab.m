function [t_enc, t_dec, C1, C2, PT] = algo_hyperchaotic_2d_matlab(data)
    % ALGO_HYPERCHAOTIC_2D - Wrapper benchmarkowy dla algorytmu 2D Chaos
    % Zgodny ze standardem testowym (algo_*.m)
    
    % 1. Przygotowanie danych
    img_orig = data.img_orig;
    img_mod  = data.img_mod;
    
    % Algorytm pracuje na double/grayscale, adapter to obsłuży.
    
    % 2. Szyfrowanie (Oryginał) - mierzymy czas i pobieramy klucze
    tic;
    [ct1, keys] = chaos_2d_adapter('encrypt', img_orig);
    t_enc = toc;
    
    % 3. Szyfrowanie (Zmodyfikowany) - tylko dla metryk NPCR/UACI
    % Klucze nas nie interesują, bo to tylko do porównania szyfrogramów
    [ct2, ~] = chaos_2d_adapter('encrypt', img_mod);
    
    % 4. Deszyfrowanie
    tic;
    pt_img = chaos_2d_adapter('decrypt', ct1, keys);
    t_dec = toc;
    
    % 5. Formatowanie wyników (wymagane uint8 dla metryk)
    C1 = uint8(ct1);
    C2 = uint8(ct2);
    PT = uint8(pt_img);
end
