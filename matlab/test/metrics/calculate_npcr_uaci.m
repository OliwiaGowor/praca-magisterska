function [NPCR, UACI] = calculate_npcr_uaci(ciphertext1, ciphertext2)
%
% This function implements the logic for NPCR and UACI measures
%
    
    % Convert to double for calculations
    ciphertext1=double(ciphertext1);
    ciphertext2=double(ciphertext2);
    
    [rows,cols]=size(ciphertext1);
    
    % Calculate NPCR (Number of Pixels Change Rate)
    NPCR=100*sum(sum(ciphertext1~=ciphertext2))/(rows*cols);
    
    % Calculate UACI (Unified Average Changing Intensity)
    UACI=100*sum(sum(abs(ciphertext1-ciphertext2)))/(rows*cols*255);

end
