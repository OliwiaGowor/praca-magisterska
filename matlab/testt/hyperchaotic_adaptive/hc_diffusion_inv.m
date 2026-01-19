function out_flat = hc_diffusion_inv(channel_flat, seq_hist)
    % HC_DIFFUSION_INV Odwraca proces dyfuzji - ZOPTIMALIZOWANA
    
    pixels = double(channel_flat);
    N = length(pixels);
    out_flat = zeros(N, 1, 'uint8');
    
    hist_counts = seq_hist(:);
    
    % Pierwszy piksel
    out_flat(1) = uint8(pixels(1));
    idx = pixels(1) + 1;
    hist_counts(idx) = hist_counts(idx) + 1;
    prev_val = pixels(1);
    
    for i = 2:N
        curr_encoded = pixels(i);
        
        % 1. Odtworzenie stanu histogramu (Szybki sort)
        [~, sorted_indices] = sort(hist_counts);
        
        % 2. Znalezienie rangi (szukamy gdzie w posortowanych jest nasza wartość)
        % curr_encoded to wartość 0..255. W sorted_indices szukamy (curr_encoded + 1).
        % Ponieważ sorted_indices to permutacja 1..256, używamy 'find' (jest szybki na małych wektorach)
        % lub można by utrzymywać odwrotną permutację, ale 'find' na 256 el. jest błyskawiczny.
        effective_rank = find(sorted_indices == (curr_encoded + 1), 1);
        
        % 3. Odtworzenie wartości z rangi (Reverse Spiral)
        if effective_rank == 1
            decoded_val = prev_val;
        else
            found_count = 1; 
            dist = 1;
            decoded_val = -1;
            
            % Pętla szukająca k-tego sąsiada
            while true
                % Lewy sąsiad
                left = prev_val - dist;
                if left >= 0
                    found_count = found_count + 1;
                    if found_count == effective_rank
                        decoded_val = left;
                        break;
                    end
                end
                
                % Prawy sąsiad
                right = prev_val + dist;
                if right <= 255
                    found_count = found_count + 1;
                    if found_count == effective_rank
                        decoded_val = right;
                        break;
                    end
                end
                dist = dist + 1;
            end
        end
        
        out_flat(i) = uint8(decoded_val);
        
        % Aktualizacja histogramu
        h_idx = curr_encoded + 1;
        hist_counts(h_idx) = hist_counts(h_idx) + 1;
        
        prev_val = curr_encoded;
    end
end