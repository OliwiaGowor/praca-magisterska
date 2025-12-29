function [NPCR, UACI] = calculate_npcr_uaci(ciphertext1, ciphertext2)
%
% This function implements the logic for NPCR and UACI measures
%

    % 8-bit representation
    C1 = uint8(ciphertext1);
    C2 = uint8(ciphertext2);

    [rows, cols] = size(C1);

    % NPCR
    NPCR = 100 * sum(C1(:) ~= C2(:)) / (rows * cols);

    % UACI
    UACI = 100 * mean(abs(double(C1(:)) - double(C2(:)))) / 255;

end

