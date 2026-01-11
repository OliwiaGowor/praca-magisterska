function [seq, params] = hc_key_gen(img)
    % HC_KEY_GEN Generuje parametry i sekwencję chaotyczną na podstawie hasha obrazu
    
    % 1. SHA-256 Hash
    try
        hashStr = DataHash(img, 'SHA-256');
    catch
        warning('Brak DataHash. Używam sumy kontrolnej jako fallback.');
        hashStr = dec2hex(sum(double(img(:)))); 
        % Dopełnienie do 64 znaków dla bezpieczeństwa logiki
        while length(hashStr) < 64
            hashStr = [hashStr, '0'];
        end
        hashStr = hashStr(1:64);
    end
    
    % Konwersja Hash (hex) na 8 bloków 32-bitowych
    blocks = uint32(zeros(1, 8));
    for i = 1:8
        hexPart = hashStr((i-1)*8 + 1 : i*8);
        blocks(i) = uint32(hex2dec(hexPart));
    end
    
    % Obliczanie parametrów (Sekcja 4.2 artykułu)
    p_alpha = bitxor(blocks(1), blocks(2));
    p_beta  = bitxor(blocks(3), blocks(4));
    p_x     = bitxor(blocks(5), blocks(6));
    p_y     = bitxor(blocks(7), blocks(8));
    
    alpha = double(p_alpha);
    beta  = double(p_beta);
    
    % Normalizacja do (0, 1]
    x_init = double(p_x) / 4294967296.0;
    y_init = double(p_y) / 4294967296.0;
    
    if x_init == 0, x_init = 0.1; end
    if y_init == 0, y_init = 0.2; end
    
    params = [alpha, beta, x_init, y_init];
    
    % Generowanie sekwencji
    [rows, cols, chans] = size(img);
    total_len = chans * (256 + rows*cols);
    
    seq = hc_map_ra2d(alpha, beta, x_init, y_init, total_len);
end