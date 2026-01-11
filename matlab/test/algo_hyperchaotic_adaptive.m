function [t_enc, t_dec, C1, C2, PT] = algo_hyperchaotic_adaptive(data)
    % ALGO_HYPERCHAOTIC_ADAPTIVE Adapter dla algorytmu 2D-RA + Adaptive Diffusion
    % Łączy framework testowy z modularną implementacją algorytmu.
    
    % Dodanie ścieżki do modułów algorytmu
    currentPath = fileparts(mfilename('fullpath'));
    algoPath = fullfile(currentPath, 'hyperchaotic_adaptive');
    addpath(algoPath);
    
    img_orig = data.img_orig;
    img_mod = data.img_mod;
    sz = data.size;
    
    %% --- SZYFROWANIE (Oryginał) ---
    tic;
    % 1. Generowanie klucza i sekwencji (zależne od hash obrazu)
    [seq_enc, params_enc] = hc_key_gen(img_orig);
    
    % 2. Właściwe szyfrowanie
    C1 = hc_encrypt(img_orig, seq_enc);
    t_enc = toc;
    
    %% --- SZYFROWANIE (Zmodyfikowany - test NPCR) ---
    % Algorytm "One-Time Pad" zależny od hasha. Dla img_mod generujemy nowy klucz.
    [seq_mod, ~] = hc_key_gen(img_mod);
    C2 = hc_encrypt(img_mod, seq_mod);
    
    %% --- DESZYFROWANIE ---
    tic;
    % Deszyfrowanie C1 przy użyciu sekwencji z img_orig
    PT = hc_decrypt(C1, seq_enc, sz);
    t_dec = toc;
    
    % Opcjonalnie: posprzątaj ścieżkę (jeśli to konieczne)
    % rmpath(algoPath); 
end