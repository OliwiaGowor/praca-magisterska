function out_flat = hc_diffusion(channel, seq_hist)
    % HC_DIFFUSION Realizuje Self-Adaptive Diffusion (Forward) - ZOPTIMALIZOWANA
    
    pixels = channel'; 
    pixels = pixels(:);
    N = length(pixels);
    out_flat = zeros(N, 1, 'uint8'); % Pre-alokacja jako uint8
    
    % Inicjalizacja histogramu sekwencją chaotyczną
    hist_counts = seq_hist(:); 
    
    % Pierwszy piksel bez zmian
    out_flat(1) = uint8(pixels(1));
    
    % Aktualizacja histogramu dla pierwszego piksela
    idx = pixels(1) + 1;
    hist_counts(idx) = hist_counts(idx) + 1;
    
    prev_val = double(pixels(1));
    
    for i = 2:N
        curr_p = double(pixels(i));
        
        % OPTYMALIZACJA 1: Zwykły sort zamiast sortrows.
        % Część ułamkowa z chaosu zapewnia unikalność, więc kolejność jest stabilna.
        % sorted_indices to bezpośrednio (wartość_piksela + 1) posortowana wg częstości.
        [~, sorted_indices] = sort(hist_counts);
        
        % 2. Obliczanie rangi wejścia (Distance logic)
        dist = abs(curr_p - prev_val);
        if curr_p == prev_val
            effective_rank = 1;
        else
            % Zoptymalizowana logika spirali (mniej if-ów)
            if curr_p < prev_val
                base_rank = 2 * dist;
                % Korekta na krawędziach
                lower = max(0, prev_val - dist + 1);
                upper = min(255, prev_val + dist - 1);
                % count_smaller = upper - lower + 1
                effective_rank = (upper - lower + 1) + 1;
            else
                base_rank = 2 * dist + 1;
                lower = max(0, prev_val - dist + 1);
                upper = min(255, prev_val + dist - 1);
                
                sibling = prev_val - dist;
                if sibling >= 0
                    effective_rank = (upper - lower + 1) + 2;
                else
                    effective_rank = (upper - lower + 1) + 1;
                end
            end
        end
        
        % 3. Mapowanie (bezpośrednio z indeksów sortowania)
        new_val_idx = sorted_indices(effective_rank);
        new_val = new_val_idx - 1; % Konwersja indeksu (1..256) na wartość (0..255)
        
        out_flat(i) = uint8(new_val);
        
        % Aktualizacja
        hist_counts(new_val_idx) = hist_counts(new_val_idx) + 1;
        prev_val = double(new_val);
    end
end