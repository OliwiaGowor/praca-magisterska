function [NPCR, UACI] = calculate_npcr_uaci(ciphertext1, ciphertext2)
%
% Ta funkcja implementuje logikę z pliku npcr_uaci_measures.m
%
    
    % Konwertuj na double dla obliczeń
    ciphertext1=double(ciphertext1);
    ciphertext2=double(ciphertext2);
    
    [rows,cols]=size(ciphertext1);
    
    % Oblicz NPCR (logika identyczna jak w pliku)
    NPCR=100*sum(sum(ciphertext1~=ciphertext2))/(rows*cols);
    
    % Oblicz UACI (logika identyczna jak w pliku)
    UACI=100*sum(sum(abs(ciphertext1-ciphertext2)))/(rows*cols*255);

end