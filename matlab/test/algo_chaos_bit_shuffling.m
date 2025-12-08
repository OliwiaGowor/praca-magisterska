function [t_enc, t_dec, C1, C2, PT] = algo_chaos_bit_shuffling(data)
    % Wrapper dla Chaos Bit Shuffling (Moysis et al.)
    
    % 1. Dane
    img_orig = data.img_orig;
    img_mod  = data.img_mod;
    sz = data.size;
    
    % Algorytm pracuje na Grayscale. Jeśli wejście jest RGB, 
    % konwersja nastąpi wewnątrz algorytmu, ale wynik będzie 2D.
    
    % 2. Szyfrowanie
    tic;
    % Wywołanie przez "dispatchera"
    [ct1, keys1] = chaos_bit_shuffling_optimized('encrypt', img_orig);
    t_enc = toc;
    
    [ct2, ~] = chaos_bit_shuffling_optimized('encrypt', img_mod);
    
    % 3. Deszyfrowanie
    tic;
    pt_gray = chaos_bit_shuffling_optimized('decrypt', ct1, keys1);
    t_dec = toc;
    
    % 4. Formatowanie
    C1 = uint8(ct1);
    C2 = uint8(ct2);
    PT = uint8(pt_gray);
end
